%% Select midfrontal components, check uniqueness and explore characteristics
% First, select the subset of significant components that are midfrontal,
% with "midfrontalness" defined by a spatial template.
% Second, perform time-frequency decomposition on components, to support 
% following analyses.
% Third, check whether the components could originate from a single neural 
% source that the GED splits out into multiple parts:
% - pairwise synchrony (should be low if unique),
% - power correlations within trials (should be low if unique),
% - power correlations across trials (should be low if unique).
% Next, explore:
% - how strongly task-modulated are components?
% - how strongly conflict-modulated are components?
% Results are saved to file for plotting script to access.
%
% Analysis code for Simon task MEEG dataset. 
%
% Authors: Marrit Zuure, Mike X Cohen
% October 2018

close all; clear;

%% Set paths
dirs = setpaths();

%% Set data import preliminaries
[sublist, ~, ~] = getICs2remove();

%% Set analysis parameters

% Midfrontal component selection threshold (higher is more template-like)
midf.r2_cutoff = 0.5;

% Time-frequency decomposition
tfd.min_freq =  2; % in Hz
tfd.max_freq = 20; % in Hz
tfd.num_frex = 40;
tfd.times2save = -700:20:1000; % in ms. Note: ITIs are 700-1200 ms.
tfd.baseline = [-500 -100]; % in ms
tfd.theta = 6; % in Hz

% Checking for uniqueness: pairwise component synchrony
compsynch.window = [-500 900]; % time window (in ms) over which to calculate synchrony (constrained by times2save and edge effects/previous trial effects)

% Checking for uniqueness: within-trial and across-trial theta amplitude correlation
thetacorr.window = [-500 900]; % time window (in ms) over which to calculate correlation (constrained by times2save and edge effects/previous trial effects)
thetacorr.alpha = 0.05 / length(sublist); % for testing temporal stability (Bonferroni-corrected)

% Task and conflict modulation
taskmod.window = [0 800]; % time window (in ms) over which to calculate modulation
conflictmod.window = [0 800]; % time window (in ms) over which to calculate modulation

%% Loop over subjects
for subno = 1:length(sublist)
    disp(['Processing subject ' num2str(subno) ' of ' num2str(length(sublist)) ' (name: ' sublist{subno} ')']);
    
    %% Load GED data
    GED_filename = [dirs.results sublist{subno} '_GED.mat'];
    load(GED_filename);
    
    %% Construct midfrontal theta template: Gaussian centered on FCz
    % Needs to be inside subject loop, because subject EEG.chanlocs
    % sometimes vary.
    
    fczidx = strcmpi('fcz',{EEG.chanlocs.labels});
    eucdist = zeros(1,EEG.nbchan);
    
    for chani = 1:EEG.nbchan
        eucdist(chani) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(fczidx).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(fczidx).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(fczidx).Z)^2 );
    end
    
    midf.template = exp(-(eucdist.^2)/(2*50^2) );
    
    %% For each significant component, construct spatial forward filter model
    %  and calculate shared variance (r^2) with template
    
    % Create EEG topoplots + compute shared variance
    midf.template_r2 = zeros(1, GED.num_comps);
    midf.ffm_EEG = zeros(EEG.nbchan, GED.num_comps);
    for c = 1:GED.num_comps
        topo = GED.evecs(1:EEG.nbchan,c)' * GED.covS(1:EEG.nbchan, 1:EEG.nbchan);
        
        % Flip component by sign of correlation with FCz template
        midf.ffm_EEG(:,c) = topo * sign(corr(topo', midf.template')); 
        
        midf.template_r2(c) = corr(topo', midf.template')^2;
    end
    
    % Create MEG topoplots
    midf.ffm_MEG = zeros(MEG.nbchan, GED.num_comps);
    
    for c = 1:GED.num_comps
        topo = GED.evecs(EEG.nbchan+1:end,c)' * GED.covS(EEG.nbchan+1:end, EEG.nbchan+1:end);
        
        % Flip MEG component by sign of max right lateral electrode
        maxval = max(abs([topo(156:186), topo(228:257)]));
        maxidx = find(abs(topo) == maxval);
        assert(length(maxidx) == 1); %checking whether there is just one electrode with that exact maximum value
        midf.ffm_MEG(:, c) = topo * sign(topo(maxidx));
    end

    %% Select midfrontal components to use in following analyses   
    comps = 1:GED.num_comps;
    midf.comps2use = comps(midf.template_r2 > midf.r2_cutoff);
    
    % For subject S10: remove first component. Later analysis finds a
    % negative task modulation, and the TF decomposition suggests that it's
    % not very theta-specific.
    if strcmpi(sublist{subno}, 'S10')
        midf.comps2use = midf.comps2use(2:end);
    end
    
    midf.num_comps = length(midf.comps2use);
    
    %% Extract just midfrontal components from component time series
    midf.compts = GED.compts(midf.comps2use,:,:);
    
    %% Time-frequency decompose midfrontal component time series
    tfd.frex = logspace(log10(tfd.min_freq),log10(tfd.max_freq),tfd.num_frex);
    tfd.tf = tfdecomp(midf.compts, MEEG, tfd.min_freq, tfd.max_freq, tfd.num_frex, 'means', tfd.times2save, tfd.baseline, trialtype);
    tfd.trials = tfdecomp(midf.compts, MEEG, tfd.theta, tfd.theta, 1, 'trials', tfd.times2save, tfd.baseline);
    
    %% Determine midfrontal component uniqueness: calculate pairwise synchrony
    % Rationale: if two components are driven by the same neural source,
    % their time series are expected to be highly synchronous.
    
    % Extract phase angle at 6 Hz from midfrontal component time series
    phase = tfdecomp(midf.compts, MEEG, tfd.theta, tfd.theta, 1, 'phase');
    
    % Time window over which to calculate synchrony
    widx = dsearchn(MEEG.times(:), compsynch.window');
    
    % List all possible pairs
    compsynch.pairs = nchoosek(1:midf.num_comps,2);
    
    % Compute vector length of mean component phase angle difference over trials
    compsynch.synch = nan(size(compsynch.pairs,1), length(MEEG.times(widx(1):widx(2))));
    for pairi = 1:size(compsynch.pairs,1)
        temp = abs(mean(exp( 1i * (phase(compsynch.pairs(pairi, 1),1,:,:) - phase(compsynch.pairs(pairi,2),1,:,:))),4));
        compsynch.synch(pairi,:) = temp(widx(1):widx(2));
    end
    
    % Calculate matrix containing mean pairwise synchrony over window
    compsynch.mtx = ones(midf.num_comps);
    for pairi = 1:size(compsynch.pairs,1)
        compsynch.mtx(compsynch.pairs(pairi,1), compsynch.pairs(pairi,2)) = mean(compsynch.synch(pairi));
        compsynch.mtx(compsynch.pairs(pairi,2), compsynch.pairs(pairi,1)) = mean(compsynch.synch(pairi));
    end
    
    %% Determine midfrontal component uniqueness: calculate within-trial power correlations
    % Rationale: if two components are driven by the same neural source,
    % their power (amplitude^2) over time should be highly correlated.
    
    % Time window over which to calculate correlations
    widx = dsearchn(tfd.times2save(:), thetacorr.window');
    
    % Calculate pairwise matrix of component amplitude correlations per trial
    thetacorr.within_trials = zeros(midf.num_comps, midf.num_comps);
    for triali = 1:MEEG.trials
        thetacorr.within_trials = thetacorr.within_trials + corr(squeeze(tfd.trials(:, 1, widx(1):widx(2), triali))');
    end
    % Normalize by number of trials
    thetacorr.within_trials = thetacorr.within_trials ./ MEEG.trials;
    
    %% Determine midfrontal component uniqueness: calculate cross-trial power correlations
    % Rationale: if two components are driven by the same neural source,
    % their mean power on a given trial should be highly correlated.
    
    % Time window over which to extract mean power
    widx = dsearchn(tfd.times2save(:), thetacorr.window');
    
    thetacorr.trial_power = squeeze(mean(tfd.trials(:,1,widx(1):widx(2),:),3));
    
    thetacorr.across_trials = corr(thetacorr.trial_power');
    
    %% Test for component temporal stability
    % Test whether component amplitudes are significantly negatively correlated:
    [~, thetacorr.p_stability] = ...
        ttest(thetacorr.across_trials(find(triu(thetacorr.across_trials,1))), 0, 'alpha', thetacorr.alpha, 'tail', 'left');
    
    %% Exploration: How strongly task-modulated are the (remaining) components?
    % Extract baseline-normalized (but not dB-scaled) theta power per trial over prestim and poststim window
    widx = dsearchn(tfd.times2save', taskmod.window');
    taskmod.score = squeeze(mean(mean(tfd.trials(:,:,widx(1):widx(2),:),3),4));
    
    %% Exploration: how strongly conflict-modulated are the (remaining) components?
    % Extract baseline-normalized (but not dB-scaled) theta power per trial
    % for cI and cC trials; compute conflict modulation as cI - cC
    
    widx = dsearchn(tfd.times2save', conflictmod.window');
    bidx = dsearchn(tfd.times2save', tfd.baseline');
    
    % Normalize cC and cI means to their own theta power baseline
    cI = mean(mean(tfd.trials(:,:,widx(1):widx(2),trialtype==2),3) - mean(tfd.trials(:,:,bidx(1):bidx(2),trialtype==2),3),4);
    cC = mean(mean(tfd.trials(:,:,widx(1):widx(2),trialtype==1),3) - mean(tfd.trials(:,:,bidx(1):bidx(2),trialtype==1),3),4);
    
    conflictmod.score = cI - cC;

    %% Save analysis results to file for later plotting and reporting
    % ...for later across-subjects analyses and plotting 
    disp('Saving results to file...');
    ana_filename = [dirs.results sublist{subno} '_ana.mat'];
    if exist(ana_filename, 'file')
        save(ana_filename, 'midf', 'tfd', 'compsynch', 'thetacorr', 'taskmod', 'conflictmod', '-append');
    else
        save(ana_filename, 'midf', 'tfd', 'compsynch', 'thetacorr', 'taskmod', 'conflictmod');
    end
    
end
disp('Run completed successfully.');