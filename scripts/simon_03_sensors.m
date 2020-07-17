%% Do sensor-level analyses + calculate sensor-level task and conflict modulation
% These analyses are split out from the rest because they require loading
% the (large) full MEEG dataset.
%
% Analysis code for Simon task MEEG dataset.
% Authors: Marrit Zuure, Mike X Cohen
% January 2019

close all; clear;

%% Set paths
dirs = setpaths();

%% Set data import/cleaning preliminaries
[sublist, EEGICs2remove, MEGICs2remove] = getICs2remove();

%% Set analysis parameters

% Sensor-level analysis task window
sensor.window = [0 800]; % in ms

% Adjacent EEG electrode synchrony
elecsynch.electrodes = {'pz', 'poz'};
% Pz and POz chosen because they're adjacent, but spatially
% distant from the expected task modulation in this dataset.

%% Loop over subjects
for subno = 1:length(sublist)
    disp(['Processing subject ' num2str(subno) ' of ' num2str(length(sublist)) '...']);
    
    %% Load GED data
    GED_filename = [dirs.results sublist{subno} '_GED.mat'];
    load(GED_filename);
   
    %% Load MEEG data
    [~, ~, MEEG.data, ~, ~] = loadMEEG(sublist, subno, EEGICs2remove, MEGICs2remove, dirs.data, 'data');
    
    %% Load previous analysis results
    ana_filename = [dirs.results sublist{subno} '_ana.mat'];
    load(ana_filename);
    
    %% Sensor-level analyses: TF decompose sensors to extract theta power in window
    sensor.tf = tfdecomp(MEEG.data, MEEG, tfd.theta, tfd.theta, 1, 'means', tfd.times2save, tfd.baseline, trialtype);
    
    %% Sensor-level analyses: extract conflict and task effect (cI - cC) at each sensor
    widx = dsearchn(tfd.times2save(:), sensor.window');
    bidx = dsearchn(tfd.times2save(:), tfd.baseline');
    
    % Normalize cI and cC per condition
    cI = mean(squeeze(sensor.tf(:,1,widx(1):widx(2),3)) - mean(sensor.tf(:,1,bidx(1):bidx(2),3),3),2);
    cC = mean(squeeze(sensor.tf(:,1,widx(1):widx(2),2)) - mean(sensor.tf(:,1,bidx(1):bidx(2),2),3),2);
    sensor.conflict = cI - cC;

    sensor.task = mean(squeeze(sensor.tf(:,1,widx(1):widx(2),1)),2);
    
    %% Calculate synchrony between adjacent electrodes
    % To use as reference value for synchrony between components.
    % Demonstrates how strong theta synchrony is that results primarily from
    % volume conduction (and thus how weak theta synchrony between components
    % actually is in comparison).
    elecsynch.window = compsynch.window;
    
    % Find electrode indices
    elec1_idx = find(strcmpi(elecsynch.electrodes{1},{MEEG.chanlocs.labels}));
    elec2_idx = find(strcmpi(elecsynch.electrodes{2},{MEEG.chanlocs.labels}));
    
    % Extract the phase angle at 6 Hz from the electrode time series
    phase = tfdecomp(MEEG.data([elec1_idx elec2_idx], :, :, :), MEEG, tfd.theta, tfd.theta, 1, 'phase');
    
    % Time window over which to calculate synchrony
    widx = dsearchn(MEEG.times(:), elecsynch.window');
    
    % Compute vector length of mean component phase angle difference over trials
    temp = abs(mean(exp( 1i * (phase(1,1,:,:) - phase(2,1,:,:))),4));
    elecsynch.synch = squeeze(temp(widx(1):widx(2)));
    
    %% Time-frequency decompose sensor FCz at all frequencies, not just theta
    fczidx = strcmpi('fcz', {MEEG.chanlocs.labels});
    sensor.fcz = tfdecomp(MEEG.data(fczidx,:,:), MEEG, tfd.min_freq, tfd.max_freq, tfd.num_frex, 'means', tfd.times2save, tfd.baseline, trialtype);
    
    %% Save results to file
    % Specifically, append to analysis results file
    disp('Saving results to file...');
    ana_filename = [dirs.results sublist{subno} '_ana.mat'];
    save(ana_filename, 'tfd', 'sensor', 'elecsynch', '-append');
end