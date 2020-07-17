%% Apply a PCA to component time courses and to component driving and
% receiving force over time, to get a sense of salient component behaviors.
%
% Analysis code for Simon task MEEG dataset.
% Author: Marrit Zuure
% January 2019

close all; clear;

%% Set paths
dirs = setpaths();

%% Set data import preliminaries
[sublist, ~, ~] = getICs2remove();

%% Set analysis parameters

% Time window for GC driver/receiver and theta power time series to go into
% PCA
pca_window = [-500 1000]; % in ms
% Broader than usual; want to include entire trial, minus edges and
% previous/next trial stimulus effects. 
% Note that component theta time series is constrained by tfd.times2save
% (which is -700:1000 ms by default).

% Frequency range for which to extract component time series TF
% decomposition
lo = 5; % in Hz
hi = 7; % in Hz

% Number of principal components for which to extract loadings
driver_PCs = 2; % based on scree plot
receiver_PCs = 2; % based on scree plot
compts_PCs = 2; % based on scree plot
MEG_PCs = 2;% based on scree plot

%% Loop over subjects
for subno = 1:length(sublist)
    disp(['Processing subject ' num2str(subno) ' of ' num2str(length(sublist)) '...']);
    
    %% Load GED data and previous analysis results
    GED_filename = [dirs.results sublist{subno} '_GED.mat'];
    load(GED_filename);
    
    ana_filename = [dirs.results sublist{subno} '_ana.mat'];
    load(ana_filename);

    %% Create data struct containing all subjects' data
     temp_struct =  struct('EEG', EEG, 'MEG', MEG, 'MEEG', MEEG, 'GED', GED, ...
        'tfd', tfd, 'thetacorr', thetacorr, 'compsynch', compsynch, 'trialtype', trialtype, ...
        'midf', midf, 'allrtidx', allrtidx, 'taskmod', taskmod, 'conflictmod', conflictmod, ...
        'gc', gc);
    
    data(subno) = temp_struct;
    clear temp_struct;
end

%% Extract variables from data struct
% Extract GC driver and receiver time courses, theta-filtered component time 
% courses, GC offsets, relative eigenvalue (serving as component time 
% series offset), task modulation, conflict modulation, midfrontalness, 
% sign of most extreme Granger interaction

gcwidx = dsearchn(gc.winmid(:), pca_window');
tswidx = dsearchn(tfd.times2save(:), pca_window');
frexidx = dsearchn(tfd.frex(:), [lo hi]');

driver_pooled = [];
receiver_pooled = [];
compts_pooled = [];
driver_baseline_pooled = [];
receiver_baseline_pooled = [];
rel_evals_pooled = [];
taskmod_pooled = [];
conflictmod_pooled = [];
driver_conflictmod_pooled = [];
receiver_conflictmod_pooled = [];
driver_taskmod_pooled = [];
receiver_taskmod_pooled = [];

for subno = 1:length(sublist)
    % Collect (theta, driver, receiver) time series per component, for PCA; de-mean
    driver_pooled = [driver_pooled; data(subno).gc.F_driver(:,gcwidx(1):gcwidx(2)) - mean(data(subno).gc.F_driver(:,gcwidx(1):gcwidx(2)),2)];
    receiver_pooled = [receiver_pooled; data(subno).gc.F_receiver(:,gcwidx(1):gcwidx(2)) - mean(data(subno).gc.F_receiver(:,gcwidx(1):gcwidx(2)),2)];
    compts_pooled = [compts_pooled; bsxfun(@minus, data(subno).tfd.tf(:,frexidx(1):frexidx(2),tswidx(1):tswidx(2),1), ...
        mean(data(subno).tfd.tf(:,frexidx(1):frexidx(2),tswidx(1):tswidx(2),1),3))];

    % Collect single values per component, for factor analysis
    driver_baseline_pooled = [driver_baseline_pooled; data(subno).gc.driver_baseline];
    receiver_baseline_pooled = [receiver_baseline_pooled; data(subno).gc.receiver_baseline];
    rel_evals_total = data(subno).GED.evals ./ sum(data(subno).GED.evals);
    rel_evals_pooled = [rel_evals_pooled; rel_evals_total(data(subno).midf.comps2use)];
    taskmod_pooled = [taskmod_pooled; data(subno).taskmod.score];
    conflictmod_pooled = [conflictmod_pooled; data(subno).conflictmod.score];
    
    widx = dsearchn(data(subno).gc.winmid(:), data(subno).gc.window');
    driver_conflictmod_pooled = [driver_conflictmod_pooled; data(subno).gc.driver_conflict];
    receiver_conflictmod_pooled = [receiver_conflictmod_pooled; data(subno).gc.receiver_conflict];
    driver_taskmod_pooled = [driver_taskmod_pooled; data(subno).gc.driver_task];
    receiver_taskmod_pooled = [receiver_taskmod_pooled; data(subno).gc.receiver_task];
end

%% Gather MEG and apply PCA

% Gather data into a matrix
ffms = {};
chanlocs = {};
for subno = 1:length(sublist)
    ffms{subno} = data(subno).midf.ffm_MEG;
    chanlocs{subno} = data(subno).MEG.chanlocs;
end
[~, new_chanlocs, new_ffms] = meanffms(ffms, chanlocs, 'labels');

% Calculate covariance matrix
MEG_cov = new_ffms * new_ffms';

% PCA
[MEG_evecs, MEG_evals] = eig(MEG_cov);

% Sort by eigenvalue
[MEG_evals, sidx] = sort(diag(MEG_evals), 'descend');
MEG_evecs = MEG_evecs(:, sidx);

%% Create covariance matrices and apply PCA
% Temporal mode; create time x time covariance matrices, not comp x comp.
% PC time series are then eigenvectors, instead of eigenvector weights
% applied to components.

% Create component theta power covariance matrix
cov_compts = zeros(size(compts_pooled,3));
% Loop over frequencies in theta range (lo to hi); take mean cov afterwards
for i = 1:size(compts_pooled,2)
    cov_compts = cov_compts + cov(squeeze(compts_pooled(:, i, :)));
end
cov_compts = cov_compts / (size(compts_pooled,2)-1);

% Apply PCA
[compts_evecs, compts_evals] = eig(cov_compts);
[compts_evals, sidx] = sort(diag(compts_evals), 'descend');
compts_evecs = compts_evecs(:,sidx);

% Create GC driver covariance matrix and apply PCA
cov_driver = driver_pooled' * driver_pooled / (size(driver_pooled,1)-1);
[driver_evecs, driver_evals] = eig(cov_driver);
[driver_evals, sidx] = sort(diag(driver_evals), 'descend');
driver_evecs = driver_evecs(:, sidx);

% Create GC receiver covariance matrix and apply PCA
cov_receiver = receiver_pooled' * receiver_pooled / (size(receiver_pooled,1)-1);
[receiver_evecs, receiver_evals] = eig(cov_receiver);
[receiver_evals, sidx] = sort(diag(receiver_evals), 'descend');
receiver_evecs = receiver_evecs(:, sidx);

%% Calculate and report component loadings on PCs
loadings = nan(size(compts_pooled,1), (compts_PCs + driver_PCs + receiver_PCs));
for i = 1:size(compts_pooled,1)
    % Calculate component loadings on different PCs: component time series
    % correlation with PC time series
    for j = 1:compts_PCs
        loadings(i,j) = corr(squeeze(mean(compts_pooled(i,:,:),2)), compts_evecs(:,j)); % mean over TF decomposition frequencies
    end
    for k = 1:driver_PCs
        loadings(i,j+k) = corr(driver_pooled(i,:)', driver_evecs(:,k));
    end
    for l = 1:receiver_PCs
        loadings(i,j+k+l) = corr(receiver_pooled(i,:)', receiver_evecs(:,l));
    end
    for m = 1:MEG_PCs
        loadings(i,j+k+l+m) = abs(corr(new_ffms(:,i), MEG_evecs(:,m))); % absolute correlation; topoplot should be sign-agnostic
    end
end

%% Flip PC eigenvectors by the most common loading sign
% If a majority of components show negative correlations with a decrease, 
% it's probably better identified as an increase

signflip = sign(sum(sign(loadings)));
% in case any variable has exactly as many positive as negative component loadings, look at the total sum of loadings
signflip(signflip==0) = sign(sum(loadings(signflip==0)));

% Flip eigenvectors
for j = 1:compts_PCs
    compts_evecs(:,j) = compts_evecs(:,j) * signflip(j);
end
for k = 1:driver_PCs
    driver_evecs(:,k) = driver_evecs(:,k) * signflip(j+k);
end
for l = 1:receiver_PCs
    receiver_evecs(:,l) = receiver_evecs(:,l) * signflip(j+k+l);
end
for m = 1:MEG_PCs
    MEG_evecs(:,m) = MEG_evecs(:,m) * signflip(j+k+l+m);
end

% Flip loadings
loadings = loadings .* signflip;

%% Save all relevant pieces of data to file
% ...for publication-quality plotting

%PCAn(alysis) because PCA is a function in MATLAB
pcan.window = pca_window;
pcan.lo = lo;
pcan.hi = hi;
pcan.driver_PCs = driver_PCs;
pcan.receiver_PCs = receiver_PCs;
pcan.compts_PCs = compts_PCs;
pcan.MEG_PCs = MEG_PCs;
pcan.driver_pooled = driver_pooled;
pcan.receiver_pooled = receiver_pooled;
pcan.compts_pooled = compts_pooled;
pcan.MEG_ffms = new_ffms;
pcan.driver_evals = driver_evals;
pcan.receiver_evals = receiver_evals;
pcan.compts_evals = compts_evals;
pcan.MEG_evals = MEG_evals;
pcan.driver_evecs = driver_evecs;
pcan.receiver_evecs = receiver_evecs;
pcan.compts_evecs = compts_evecs;
pcan.MEG_evecs = MEG_evecs;

disp('Saving results to file...');
filename = [dirs.results 'pca.mat']; % not subject-specific so can't save to subject file
save(filename, 'pcan');
