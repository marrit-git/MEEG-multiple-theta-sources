%% Check how much noise is needed in the eigenvectors to get component amplitude correlations to be as low as they are
% ...demonstrating that it's unlikely that the components are a single
% source split out into multiple parts because of noise in the
% recordings.
%
% Analysis code for Simon task MEEG dataset.
% Author: Marrit Zuure

close all; clear;

%% Set paths
dirs = setpaths();

%% Set data import/cleanup preliminaries
[sublist, EEGICs2remove, MEGICs2remove] = getICs2remove();

%% Loop over subjects
for subno = 1:length(sublist)
    disp(['Processing subject ' num2str(subno) ' of ' num2str(length(sublist)) '...']);
    
    %% Load GED and analysis data
    GED_filename = [dirs.results sublist{subno} '_GED.mat'];
    load(GED_filename);
    
    ana_filename = [dirs.results sublist{subno} '_ana.mat'];
    load(ana_filename);
    
    %% Load EEG/MEG data
    [~, ~, MEEG.data, allrtidx, trialtype] = loadMEEG(sublist, subno, EEGICs2remove, MEGICs2remove, dirs.data, 'data');
    
    %% Create noisy eigenvectors; calculate time series; correlate amplitudes
    
    iterations = 50;
    % 50 noisy children per eigenvector produces (50*50) / 2 - 50 = 1200
    % pairwise correlations with unique pairs.
    
    noise = 0:0.1:10;
    
    noisecorr_within = zeros(midf.num_comps, length(noise), iterations, iterations);
    noisecorr_across = zeros(midf.num_comps, length(noise), iterations, iterations);
    
    for noisei = 1:length(noise)
        for c = 1:midf.num_comps
            % Take parent component; create 100 noisy copies
            parent = GED.evecs(:,midf.comps2use(c)); % component eigenvector
            children = repmat(parent, [1, iterations]); 
            children = children + randn(size(children)) * std(parent) * noise(noisei); % add Gaussian noise to each child
            
            % Recalculate component time series using noisy
            % eigenvectors
            noisyts = nan(iterations, MEEG.pnts, MEEG.trials);
            for i = 1:iterations
                temp = children(:,i)' * reshape(MEEG.data, MEEG.nbchan, []);
                noisyts(i,:,:) = reshape(temp, MEEG.pnts, MEEG.trials);
            end
            
            % Time-frequency decompose at theta
            noisytfd = tfdecomp(noisyts, MEEG, 6, 6, 1, 'trials', tfd.times2save, tfd.baseline);
            
            % Calculate pairwise correlations between children -
            % within-trial
            widx = dsearchn(tfd.times2save(:), thetacorr.window');
            for triali = 1:MEEG.trials
                noisecorr_within(c,noisei,:,:) = squeeze(noisecorr_within(c,noisei,:,:)) + corr(squeeze(noisytfd(:,1,widx(1):widx(2), triali))');
            end
            % Normalize by number of trials
            noisecorr_within(c,noisei,:,:) = noisecorr_within(c,noisei,:,:) ./ MEEG.trials;
            
            % Calculate pairwise correlations between children -
            % across trials
            noisy_amplitude = squeeze(mean(noisytfd(:,1,:,:),3));
            noisecorr_across(c,noisei,:,:) = corr(noisy_amplitude');
        end
    end
    noisecorr.within = noisecorr_within;
    noisecorr.across = noisecorr_across;
    
    %% T-test against actual component correlations
    %(which is something that should be done per subject, not across subjects)
    
    hw = nan(1, length(noise));
    ha = nan(1, length(noise));
    pw = zeros(1, length(noise));
    pa = zeros(1, length(noise));
    % Test at which amount of added noise the noise r is no longer
    % significantly higher than the component r.
    % Intentionally not applying multiple comparisons correction; trying to
    % get an upper bound.
    for noisei = 1:length(noise)
        triu_idx = logical(triu(ones(iterations, iterations), 1));
        poolw = reshape(squeeze(noisecorr.within(:, noisei, triu_idx)), 1, []);
        poola = reshape(squeeze(noisecorr.across(:, noisei, triu_idx)), 1, []);
        [hw(noisei), pw(noisei), ~, wstats(noisei)] = ttest2(thetacorr.within_trials(find(triu(thetacorr.within_trials,1))), poolw, 'tail', 'left');
        [ha(noisei), pa(noisei), ~, astats(noisei)] = ttest2(thetacorr.across_trials(find(triu(thetacorr.across_trials,1))), poola, 'tail', 'left');
    end
    
    noisecorr.hw = hw;
    noisecorr.pw = pw;
    noisecorr.wstats = wstats;
    noisecorr.ha = ha;
    noisecorr.pa = pa;
    noisecorr.astats = astats;
    noisecorr.noise = noise;
    noisecorr.iterations = iterations;
    
    % For checking
    disp('Hw:');
    disp(hw);
    disp('Pw:');
    disp(pw);
    disp('Ha:');
    disp(ha)
    disp('Pa:');
    disp(pa);

    %% Save analysis results to file for later plotting and reporting
    % ...for later across-subjects analyses and plotting

    disp('Saving results to file...');
    ana_filename = [dirs.results sublist{subno} '_ana.mat'];
    save(ana_filename, 'noisecorr', '-append');
end
disp('Run completed successfully.');