%% Perform generalized eigendecomposition and determine component significance
% using permutation testing.
%
% Analysis code for Simon task MEEG dataset. 
%
% Authors: Marrit Zuure, Mike X Cohen
% October 2018

close all; clear;

%% Set paths
dirs = setpaths();

%% Set data import/cleaning preliminaries
[sublist, EEGICs2remove, MEGICs2remove] = getICs2remove();

%% Set GED and permutation test parameters

% GED filter settings
filter_lo = 5; % lower theta bound in Hz
filter_hi = 7; % upper theta bound in Hz
filter_order = 3000; % higher order is better frequency resolution, poorer temporal resolution

% GED time window (in ms) to use in constructing covariance matrices
covS_window = [0 800]; % signal matrix
covR_window = [0 800]; % noise matrix

% Permutation testing
num_perm = 1000; % 1000 recommended for accuracy, but will be SLOW (take ~1,5 days)

%% Loop over subjects
for subno = 1:length(sublist)
    disp(['Processing subject ' num2str(subno) ' of ' num2str(length(sublist)) '...']);
    
    %% Load EEG/MEG data
    [EEG, MEG, MEEG, allrtidx, trialtype] = loadMEEG(sublist, subno, EEGICs2remove, MEGICs2remove, dirs.data, 'all');

    %% Construct frequency filter kernel
    nyquist = MEEG.srate/2;
    
    try
        kernel = fir1(filter_order, [filter_lo/nyquist filter_hi/nyquist]);
    catch
        kernel = fir1(filter_order, [filter_lo/nyquist filter_hi/nyquist]);
    end
    % This fails the first time but completes the second time. Seems to happen
    % because EEGlab overrides fir1 with an Octave function, which is removed
    % from path when MATLAB detects this.
    
    %% Perform a generalized eigendecomposition on theta-filtered vs. unfiltered data
    [evals, evecs, evecs_rel, covS, covR, evecs_unnormed] = performGED(MEEG, covS_window, covR_window, kernel);
    orig_evals = evals;
    orig_evecs = evecs_unnormed;

    %% Check for existing permutation testing data; if not, run permutation test
    permtest_filename = [dirs.results sublist{subno} '_' num2str(num_perm) '_permtest.mat'];
    if exist(permtest_filename,'file')
        load(permtest_filename, 'perm_settings');
        % Check whether loaded permutation data were generated using the 
        % same settings (filter, window) as the GED. If not, the Z-scoring 
        % and thus the significance thresholding will not be valid.
        % Obviously, MEEG also needs to be the same, but we don't check for
        % that.
        assert(all(perm_settings == [covS_window, covR_window, kernel]), ...
            [ permtest_filename ' may not have been generated using the same ' ...
            'filter kernel and time windows as were used for the GED.']);
        load(permtest_filename);
    else
        [perm_evals, perm_settings, perm_evecs, perm_covS] = permuteGED(MEEG, covS_window, covR_window, kernel, num_perm, '95');
        save(permtest_filename, 'perm_evals', 'perm_settings', 'perm_evecs', 'perm_covS');
    end
    
    %% Determine significant components
    num_comps = sum(evals > perm_evals');

    %% Cleanup: Drop components that potentially represent eigenplanes
    
    % Flag significant eigenvalues that are closer together than 1% 
    % as possible repeat eigenvalues (indicating eigenplanes).
    repeat_evals_idx = zeros(1,num_comps);
    for i = 1:num_comps-1
        tolerance = evals(i) / 100;
        repeat_evals_idx(i+1) = evals(i+1) > evals(i)-tolerance;
    end
    repeat_evals_idx = logical(repeat_evals_idx);

    % Drop components belonging to repeating eigenvalues (second component
    % of the pair only)
    comps2drop = find(repeat_evals_idx);
    repeat_evals = evals(repeat_evals_idx); % replace logical index with actual eigenvalues
    for i = 1:sum(repeat_evals_idx)
        evals = [evals(1:comps2drop(i)-1); evals(comps2drop(i)+1:end)];
        evecs = [evecs(:,1:comps2drop(i)-1), evecs(:,comps2drop(i)+1:end)];
        evecs_rel = [evecs_rel(:,1:comps2drop(i)-1), evecs_rel(:,comps2drop(i)+1:end)];
    end
    num_comps = num_comps - sum(repeat_evals_idx);
    
    %% Construct time series for components above significance threshold
    compts = zeros(num_comps, MEEG.pnts, MEEG.trials);

    % Note: component time series are scaled by relative eigenvalue,
    % facilitating amplitude comparisons within and between subjects.
    for c = 1:num_comps
        compts_temp = evecs_rel(:,c)' * reshape(MEEG.data, MEEG.nbchan, []); % apply sensor weightings to concatenated trials
        compts(c,:,:) = reshape(compts_temp, MEEG.pnts, MEEG.trials); % reshape back to individual trials
    end

    %% Save GED results to file
    
    % First, trim data from EEG, MEG, MEEG structs, to prevent saving and
    % loading from taking impractically long. The other analysis files tend
    % to primarily use the compts and the other struct fields anyway. Full
    % data sets can be loaded in when needed by calling loadMEEG with the 
    % 'data' flag.
    EEG = rmfield(EEG, 'data');
    MEG = rmfield(MEG, 'data');
    MEEG = rmfield(MEEG, 'data');
    
    % Second, put results from GED in struct, to distinguish between
    % variables intrinsic to data and variables following from GED.
    GED.orig_evals = orig_evals; % for plotting
    GED.orig_evecs = orig_evecs; % for plotting
    GED.evals = evals;
    GED.repeat_evals = repeat_evals;
    GED.repeat_evals_idx = repeat_evals_idx;
    GED.evecs = evecs;
    GED.num_comps = num_comps;
    GED.compts = compts;
    GED.num_perm = num_perm;
    GED.perm_evals = perm_evals;
    GED.covS = covS; % for plotting
    GED.covR = covR; % for plotting
    GED.kernel = kernel;
    GED.lo = filter_lo;
    GED.hi = filter_hi;
    
    disp('Saving results to file...');
    GED_filename = [dirs.results sublist{subno} '_GED.mat'];
    save(GED_filename, 'EEG', 'MEG', 'MEEG', 'allrtidx', ...
        'trialtype', 'GED');
end
disp('Run completed successfully.');