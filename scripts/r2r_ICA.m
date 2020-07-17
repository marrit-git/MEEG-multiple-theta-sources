%% Check robustness of components by applying ICA to theta-filtered data.
% Analysis code for Simon task MEEG dataset, newly added during revision.
% Author: Marrit Zuure
% June 2019

close all; clear;

dirs = setpaths();

% Load EEGlab
eeglab;

%% Load MEEG data for all subjects

[sublist, EEGICs2remove, MEGICs2remove] = getICs2remove();

% Create variable to store results for all subjects
MEEGica = {};

for subno = 1:length(sublist)

    disp(['Loading data for subject ' num2str(subno) '...']);
    [EEG, MEG, MEEG, allrtidx, trialtype] = loadMEEG(sublist, subno, EEGICs2remove, MEGICs2remove, dirs.data, 'all');

    %% Filter data in theta band (5-7 Hz)

    filter_lo = 5; % lower theta bound in Hz
    filter_hi = 7; % upper theta bound in Hz
    filter_order = 3000; % higher order is better frequency resolution, poorer temporal resolution

    % Construct frequency filter kernel
    nyquist = MEEG.srate/2;

    try
        kernel = fir1(filter_order, [filter_lo/nyquist filter_hi/nyquist]);
    catch
        kernel = fir1(filter_order, [filter_lo/nyquist filter_hi/nyquist]);
    end

    disp('Filtering data...');
    fft_len = size(MEEG.data,2)+length(kernel)-1; % number of time points/frequencies for fft to return
    trim_len = (length(kernel)-1)/2; % number of time points to trim from start and end of result
    fdata = 2*real( ifft( bsxfun(@times,fft(MEEG.data,fft_len,2),fft(kernel, fft_len)) ,[],2) );
    fdata = reshape(fdata(:,trim_len+1:end-trim_len,:),size(MEEG.data,1), MEEG.pnts, MEEG.trials);

    %% Perform ICA using EEGlab

    window = [0 800]; % Time window (in ms) to use
    widx = dsearchn(MEEG.times(:), window');

    imported_EEG = EEG; % store original EEG data for safekeeping; EEGlab ICA needs input data to be named EEG
    EEG = MEEG;
    EEG.data = fdata(:, widx(1):widx(2), :); % data that goes into ICA

    % Add some values so EEGlab accepts this as a proper EEG struct
    EEG.pnts = size(EEG.data,2);
    EEG.icaweights = [];
    EEG.setname = [];
    EEG.icawinv = [];
    EEG.icasphere = [];
    EEG.icaact = [];
    EEG.xmin = window(1) / 1000; % in s
    EEG.xmax = window(2) / 1000; % in s

    % Apply ICA algorithm
    EEG = pop_runica(EEG, 'icatype', 'jader', 'dataset', 1, 'options', {60});
    % Look for 60 sources; roughly upper bound of significant components across subjects

    % Restore variable names
    MEEGica{subno}.icaweights = EEG.icaweights;
    MEEGica{subno}.icaact = EEG.icaact;
    MEEGica{subno}.icawinv = EEG.icawinv;
    MEEGica{subno}.icachansind = EEG.icachansind;

    EEG = imported_EEG;
    clear imported_EEG;

end

%% Check: >1 component (pref. near top) that are midfrontal in EEG and differ in MEG?

for subno = 1:length(sublist)
    load([dirs.results sublist{subno} '_GED.mat']); % Load file locally
    
    % Create midfrontal EEG template
    fczidx = strcmpi('fcz',{EEG.chanlocs.labels});
    eucdist = zeros(1,EEG.nbchan);
    for chani = 1:EEG.nbchan
        eucdist(chani) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(fczidx).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(fczidx).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(fczidx).Z)^2 );
    end
    midf_template = exp(-(eucdist.^2)/(2*50^2) );
    % Plot as sanity check:
    % figure; topoplot(midf_template, EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');

    % Sort ICs by shared variance with midfrontal template
    for i = 1:60
        temp = corr(MEEGica{subno}.icawinv(1:EEG.nbchan,i), midf_template');
        flip(i) = sign(temp); % Use this to flip components that are negatively correlated with midfrontal template
        R_squared(i) = temp^2;
    end
    [~, sidx] = sort(R_squared, 'descend');

    % Count the number of midfrontal ICs
    disp(['Subject ' sublist{subno} ' has ' num2str(sum(R_squared > 0.5)) ' midfrontal ICs (R^2 > 0.5)']);

    % Plot 30 ICs sorted by midfrontalness (out of top 60)
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    eegplots = [1:10 21:30 41:50];
    megplots = eegplots+10;
    for i = 1:30
        IC2plot = sidx(i);
        topo1 = subplot(6,10,eegplots(i));
        title(['EEG IC ' num2str(IC2plot)]);
        topoplot(flip(IC2plot)*MEEGica{subno}.icawinv(1:EEG.nbchan, IC2plot), EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');

        topo2 = subplot(6,10,megplots(i));
        title(['MEG IC ' num2str(IC2plot)]);
        topoplot(flip(IC2plot)*MEEGica{subno}.icawinv(EEG.nbchan+1:end, IC2plot), MEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
        axis([axis * 1.2]);

        clim1 = max(abs(caxis(topo1)));
        clim2 = max(abs(caxis(topo2)));
        clim = max(clim1, clim2);
        caxis(topo1, [-1 1]*clim);
        caxis(topo2, [-1 1]*clim);
    end
    colormap(bluewhitered);
    suptitle(['30 most midfrontal ICs (from top 60 theta ICs) for subject ' sublist{subno} '; ' num2str(sum(R_squared > 0.5)) ' R^2 > 0.5']);

end