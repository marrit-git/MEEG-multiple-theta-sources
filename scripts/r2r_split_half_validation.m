%% Perform component robustness analysis by performing GED on subsampled data
% and identifying which components remain intact.
% Analysis code for Simon task MEEG dataset, newly added during revision.
% Author: Marrit Zuure
% June 2019

close all; clear;

dirs = setpaths();

% Load EEGlab
eeglab;

[sublist, EEGICs2remove, MEGICs2remove] = getICs2remove();

%% Load data per subject
for subno = 1:length(sublist) 
    disp(['Loading data for subject ' num2str(subno) '...']);
    [EEG, MEG, MEEG, allrtidx, trialtype] = loadMEEG(sublist, subno, EEGICs2remove, MEGICs2remove, dirs.data, 'all');
    
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

    %% Perform split-half validation a la Groppe, Makeig, & Kutas (2009). 
    % First, split the data in equal halves, alternating trials

    data1 = MEEG.data(:,:,1:2:end);
    data2 = MEEG.data(:,:,2:2:end);

    fdata1 = fdata(:,:,1:2:end);
    fdata2 = fdata(:,:,2:2:end);

    %% Compute covariance matrices per trial

    window = [0 800]; % Time window (in ms) to use
    widx = dsearchn(MEEG.times(:), window');

    % First half
    covS1 = zeros(MEEG.nbchan);
    for i = 1:size(data1,3)
        temp = zscore(squeeze(fdata1(:, widx(1):widx(2), i)));
        covS1 = covS1 + temp*temp' / MEEG.nbchan;
    end
    covS1 = covS1 / i;

    covR1 = zeros(MEEG.nbchan);
    for i = 1:size(data1,3)
        temp = zscore(squeeze(data1(:, widx(1):widx(2), i)));
        covR1 = covR1 + temp*temp' / MEEG.nbchan;
    end
    covR1 = covR1 / i;

    % Regularization
    g = 0.01;
    covRr1 = (1-g)*covR1 + g*mean(eig(covR1))*eye(MEEG.nbchan);

    % Second half
    covS2 = zeros(MEEG.nbchan);
    for i = 1:size(data2,3)
        temp = zscore(squeeze(fdata2(:, widx(1):widx(2), i)));
        covS2 = covS2 + temp*temp' / MEEG.nbchan;
    end
    covS2 = covS2 / i;

    covR2 = zeros(MEEG.nbchan);
    for i = 1:size(data2,3)
        temp = zscore(squeeze(data2(:, widx(1):widx(2), i)));
        covR2 = covR2 + temp*temp' / MEEG.nbchan;
    end
    covR2 = covR2 / i;

    % Regularization
    g = 0.01;
    covRr2 = (1-g)*covR2 + g*mean(eig(covR2))*eye(MEEG.nbchan);

    %% Save covariance matrices to file (server)

    save([dirs.results 'r2r_' sublist{subno} '_split_half.mat'], 'covS1', 'covR1', 'covRr1', 'covS2', 'covR2', 'covRr2');
end

%% Perform GEDs on both data halves
% Create struct array to store variables
splithalf = struct('pairedcomps',{},'eegtopos',{},'megtopos',{},'num_comps',{});

for subno = 1:length(sublist)
    disp(['Computing GEDs and null distribution for subject ' num2str(subno) '...']);
    
    % Load covariance matrices (local)
    load([dirs.results 'r2r_' sublist{subno} '_split_half.mat']);
    load([dirs.results sublist{subno} '_GED.mat']);

    % Perform GEDs
    [evecs1, evals1] = eig(covS1, covRr1);
    [evals1, sidx] = sort(diag(evals1), 'descend');
    evecs1 = evecs1(:, sidx);

    [evecs2, evals2] = eig(covS2, covRr2);
    [evals2, sidx] = sort(diag(evals2), 'descend');
    evecs2 = evecs2(:, sidx);

    % Norm eigenvectors to unit length
    for v = 1:size(evecs1,2)
        evecs1(:,v) = evecs1(:,v)/norm(evecs1(:,v));
        evecs2(:,v) = evecs2(:,v)/norm(evecs2(:,v));
        GED.orig_evecs(:,v) = GED.orig_evecs(:,v)/norm(GED.orig_evecs(:,v)); % these vectors weren't normalized yet
    end

    %% Compute topographies for full set and both halves
    % Limiting ourselves to top 60 components (= roughly the latest significant component among subjects + 10)

    eegtopos = [];
    megtopos = [];
    eegtopos1 = [];
    megtopos1 = [];
    eegtopos2 = [];
    megtopos2 = [];
    for i = 1:60
        eegtopos(i,:) = GED.orig_evecs(1:EEG.nbchan,i)' * GED.covS(1:EEG.nbchan, 1:EEG.nbchan);
        megtopos(i,:) = GED.orig_evecs(EEG.nbchan+1:end,i)' * GED.covS(EEG.nbchan+1:end, EEG.nbchan+1:end);
        eegtopos1(i,:) = evecs1(1:EEG.nbchan,i)' * covS1(1:EEG.nbchan, 1:EEG.nbchan);
        megtopos1(i,:) = evecs1(EEG.nbchan+1:end,i)' * covS1(EEG.nbchan+1:end, EEG.nbchan+1:end);
        eegtopos2(i,:) = evecs2(1:EEG.nbchan,i)' * covS2(1:EEG.nbchan, 1:EEG.nbchan);
        megtopos2(i,:) = evecs2(EEG.nbchan+1:end,i)' * covS2(EEG.nbchan+1:end, EEG.nbchan+1:end);
    end

    %% Create null distribution
    topodist = @(eeg1, meg1, eeg2, meg2) ...
        1 - ...
        abs( ... % absolute
        mean( ... % first part to take mean over ( = EEG)
        [(eeg1 * eeg2') ... % EEG dot product
        / (norm(eeg1)*norm(eeg2)) ... % EEG topo magnitudes
        , ... % second part to take mean over (= MEG)
        (meg1 * meg2') ...
        / (norm(meg1)*norm(meg2))] ... % MEG topo magnitudes
        ));
    
    all_eegtopos = [];
    all_megtopos = [];
    all_eegtopos1 = [];
    all_megtopos1 = [];
    all_eegtopos2 = [];
    all_megtopos2 = [];
    for i = 1:MEEG.nbchan
        all_eegtopos(i,:) = GED.orig_evecs(1:EEG.nbchan,i)' * GED.covS(1:EEG.nbchan, 1:EEG.nbchan);
        all_megtopos(i,:) = GED.orig_evecs(EEG.nbchan+1:end,i)' * GED.covS(EEG.nbchan+1:end, EEG.nbchan+1:end);
        all_eegtopos1(i,:) = evecs1(1:EEG.nbchan,i)' * covS1(1:EEG.nbchan, 1:EEG.nbchan);
        all_megtopos1(i,:) = evecs1(EEG.nbchan+1:end,i)' * covS1(EEG.nbchan+1:end, EEG.nbchan+1:end);
        all_eegtopos2(i,:) = evecs2(1:EEG.nbchan,i)' * covS2(1:EEG.nbchan, 1:EEG.nbchan);
        all_megtopos2(i,:) = evecs2(EEG.nbchan+1:end,i)' * covS2(EEG.nbchan+1:end, EEG.nbchan+1:end);
    end
    
    nulldistr = nan(MEEG.nbchan, MEEG.nbchan, 3);
    for i = 1:MEEG.nbchan
        for j = 1:MEEG.nbchan
            nulldistr(i,j,1) = topodist(all_eegtopos(i,:), all_megtopos(i,:), all_eegtopos1(j,:), all_megtopos1(j,:));
            nulldistr(i,j,2) = topodist(all_eegtopos(i,:), all_megtopos(i,:), all_eegtopos2(j,:), all_megtopos2(j,:));
            nulldistr(i,j,3) = topodist(all_eegtopos1(i,:), all_megtopos1(i,:), all_eegtopos2(j,:), all_megtopos2(j,:));
        end
    end
    
    nulldistr_file = [dirs.results 'r2r_splithalf_nulldistr.mat'];
    if exist(nulldistr_file, 'file')
        load(nulldistr_file);
        nulldistr_all = [nulldistr_all, nulldistr(:)'];
        save(nulldistr_file, 'nulldistr_all');
    else
        nulldistr_all = nulldistr(:)';
        save(nulldistr_file, 'nulldistr_all');
    end

    %% Create midfrontal EEG template
    fczidx = strcmpi('fcz',{EEG.chanlocs.labels});
    eucdist = zeros(1,EEG.nbchan);
    for chani = 1:EEG.nbchan
        eucdist(chani) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(fczidx).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(fczidx).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(fczidx).Z)^2 );
    end
    midf_template = exp(-(eucdist.^2)/(2*50^2) );

    midfrontalness = nan(60,3);
    flip = nan(60,3);
    for i = 1:60
        temp = corr(eegtopos(i,1:EEG.nbchan)', midf_template');
        midfrontalness(i,1) = temp^2;
        flip(i,1) = sign(temp);
        temp = corr(eegtopos1(i,1:EEG.nbchan)', midf_template');
        midfrontalness(i,2) = temp^2;
        flip(i,2) = sign(temp);
        temp = corr(eegtopos2(i,1:EEG.nbchan)', midf_template');
        midfrontalness(i,3) = temp^2;
        flip(i,3) = sign(temp);
    end

    eegtopos = eegtopos .* flip(:,1);
    eegtopos1 = eegtopos1 .* flip(:,2);
    eegtopos2 = eegtopos2 .* flip(:,3);
    megtopos = megtopos .* flip(:,1);
    megtopos1 = megtopos1 .* flip(:,2);
    megtopos2 = megtopos2 .* flip(:,3);

    [midf, sidx] = sort(midfrontalness(:,1), 'descend');
    [~, sidx1] = sort(midfrontalness(:,2), 'descend');
    [~, sidx2] = sort(midfrontalness(:,3), 'descend');

    % % Sanity check: plot top 25 EEG topos for each, sorted by midfrontalness
    % figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    % for i = 1:25
    %     subplot(5,5,i);
    %     topoplot(eegtopos(sidx(i),:), EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
    % end
    % suptitle('topos full set');
    % colormap bluewhitered;
    % 
    % saveas(gcf, [dirs.plots 'r2r_03_fullset.png']);
    % close;
    % 
    % figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    % for i = 1:25
    %     subplot(5,5,i);
    %     topoplot(eegtopos1(sidx1(i),:), EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
    % end
    % suptitle('topos split half 1');
    % colormap bluewhitered;
    % 
    % saveas(gcf, [dirs.plots 'r2r_03_splithalf1.png']);
    % close;
    % 
    % figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    % for i = 1:25
    %     subplot(5,5,i);
    %     topoplot(eegtopos2(sidx2(i),:), EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
    % end
    % suptitle('topos split half 2');
    % colormap bluewhitered;
    % 
    % saveas(gcf, [dirs.plots 'r2r_03_splithalf2.png']);
    % close;

    %% Compute topographical distance between full set and both halves

    % New distance definition: compute dot products for EEG and MEG separately
    % (otherwise MEG is weighted more in the distance metric). Take mean.
    topodist = @(eeg1, meg1, eeg2, meg2) ...
        1 - ...
        abs( ... % absolute
        mean( ... % first part to take mean over ( = EEG)
        [(eeg1 * eeg2') ... % EEG dot product
        / (norm(eeg1)*norm(eeg2)) ... % EEG topo magnitudes
        , ... % second part to take mean over (= MEG)
        (meg1 * meg2') ...
        / (norm(meg1)*norm(meg2))] ... % MEG topo magnitudes
        ));

    dist1 = nan(60,60);
    dist2 = nan(60,60);
    for i = 1:60
        for j = 1:60
            dist1(i,j) = topodist(eegtopos(i,:), megtopos(i,:), eegtopos1(j,:), megtopos1(j,:));
            dist2(i,j) = topodist(eegtopos(i,:), megtopos(i,:), eegtopos2(j,:), megtopos2(j,:));
        end
    end

    % Make backup of variables; these are changed in a later loop
    full_dist1 = dist1;
    full_dist2 = dist2;

    %% Greedily pair together triplets by topography

    % Start from full set
    pairedcomps = nan(60,5);
    dist1 = full_dist1;
    dist2 = full_dist2;
    for i = 1:60
        topos1_dist = min(dist1(:));
        [pairedcomps(i,1), pairedcomps(i,2)] = find(ismember(dist1, min(dist1(:)))); % get most similar: topos + topos1
        topos2_dist = min(dist2(pairedcomps(i,1),:));
        [temp, pairedcomps(i,3)] = find(ismember(dist2, min(dist2(pairedcomps(i,1),:)))); % get most similar: current component in topos + topos2
        assert(temp == pairedcomps(i,1));
        pairedcomps(i,4) = topos1_dist;
        pairedcomps(i,5) = topos2_dist;

        % NaN these components in the dist matrices (leaving row and column
        % indexing intact)
        dist1(pairedcomps(i,1), :) = nan;
        dist1(:, pairedcomps(i,2)) = nan;
        dist2(pairedcomps(i,1), :) = nan;
        dist2(:, pairedcomps(i,3)) = nan;
    end

    % Assert that components in each column are unique (no double pairing)
    for i = 1:3
        assert(length(pairedcomps(:,i)) == length(unique(pairedcomps(:,i))));
    end
    
    % Save variables in struct for next loop
    splithalf{subno}.pairedcomps = pairedcomps;
    splithalf{subno}.eegtopos = eegtopos;
    splithalf{subno}.megtopos = megtopos;
    splithalf{subno}.midf_template = midf_template;
    splithalf{subno}.num_comps = GED.num_comps;
    splithalf{subno}.EEG_chanlocs = EEG.chanlocs;
    splithalf{subno}.MEG_chanlocs = MEG.chanlocs;
end

for subno = 1:length(sublist)
    
    % Extract variables from struct
    pairedcomps = splithalf{subno}.pairedcomps;
    eegtopos = splithalf{subno}.eegtopos;
    megtopos = splithalf{subno}.megtopos;
    num_comps = splithalf{subno}.num_comps;
    EEG_chanlocs = splithalf{subno}.EEG_chanlocs;
    MEG_chanlocs = splithalf{subno}.MEG_chanlocs;
    midf_template = splithalf{subno}.midf_template;
    
    %% Approximate null hypothesis distributions for the distances by computing the 3n2 all-to-all component distances for all sets of components
    nulldistr_file = [dirs.results 'r2r_splithalf_nulldistr.mat'];
    load(nulldistr_file);
    nulldistr = nulldistr_all;

%     % Plot as sanity check
%     figure('units', 'normalized', 'outerposition', [0 0 1 1]);
%     hist(nulldistr(:), 30);
%     title('Null distribution of topographical distances between all 328*328*3 split-half components');
%     saveas(gcf, [dirs.plots 'r2r_nulldist.png']);
%     close;

    % Get topographies with 95th percentile similarity ( = 5th percentile distance)
    cutoff = prctile(nulldistr(:), 5);

    % Select the previously paired triplets for which none of the distances
    % exceed the cutoff
    idx = pairedcomps(:,4) < cutoff & pairedcomps(:,5) < cutoff;
    pairedcomps = pairedcomps(idx,:);
    fprintf('Rejected %i out of 60 triplets for exceeding 5th percentile distance cutoff of %.3f\n', 60 - length(pairedcomps), cutoff);

    % Remove the previously paired triplets for which the full set component is
    % outside the significance range
    idx = pairedcomps(:,1) < num_comps;
    pairedcomps = pairedcomps(idx,:);
    fprintf('Rejected %i out of %i triplets for not being significant in the full set (i > %i).', length(idx) - length(pairedcomps), length(idx), GED.num_comps);

    %% Plot homologous topographies/report on homologous vs. non-homologous components

    % Sort remaining (paired) components by midfrontalness; plot most
    % midfrontal
    midfrontalness = zeros(size(pairedcomps,1), 1);
    for i = 1:size(pairedcomps,1)
        midfrontalness(i) = corr(eegtopos(pairedcomps(i,1),:)', midf_template')^2;
    end
    [midfrontalness_sorted, sidx] = sort(midfrontalness, 'descend');

    % Define function to compute p-value from null distribution + given value
    pvalue = @(nulldist, value) (sum(nulldist < value) + 0.5*sum(nulldist == value)) / length(nulldist);

    eeglocs = [1:2:10];
    meglocs = [2:2:10];
    eeglocs1 = eeglocs + 10;
    meglocs1 = meglocs + 10;
    eeglocs2 = eeglocs + 20;
    meglocs2 = meglocs + 20;

    p = 0;
    for j = 1:size(pairedcomps,1)
        i = mod(j-1,5) + 1;
        if i == 1
            figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        end

        topo1 = subplot(3, 10, eeglocs(i));
        topoplot(eegtopos(pairedcomps(sidx(j),1),:), EEG_chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp'); % full set, EEG
        title(['Comp ' num2str(pairedcomps(sidx(j),1)) ', R^2 = ' num2str(midfrontalness_sorted(j), '%.4f')])
        topo2 = subplot(3, 10, meglocs(i));
        topoplot(megtopos(pairedcomps(sidx(j),1),:), MEG_chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp'); % full set, MEG
        axis([axis * 1.2]);
        clim1 = max(abs(caxis(topo1)));
        clim2 = max(abs(caxis(topo2)));
        clim = max(clim1, clim2);
        caxis(topo1, [-1 1]*clim);
        caxis(topo2, [-1 1]*clim);

        topo1 = subplot(3, 10, eeglocs1(i));
        topoplot(eegtopos1(pairedcomps(sidx(j),2),:), EEG_chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp'); % full set, EEG
        title(['Comp ' num2str(pairedcomps(sidx(j),2)) ', p = ' num2str(pvalue(nulldistr(:), pairedcomps(sidx(j),4)), '%.3f')])
        topo2 = subplot(3, 10, meglocs1(i));
        topoplot(megtopos1(pairedcomps(sidx(j),2),:), MEG_chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp'); % full set, MEG
        axis([axis * 1.2]);
        clim1 = max(abs(caxis(topo1)));
        clim2 = max(abs(caxis(topo2)));
        clim = max(clim1, clim2);
        caxis(topo1, [-1 1]*clim);
        caxis(topo2, [-1 1]*clim);

        topo1 = subplot(3, 10, eeglocs2(i));
        topoplot(eegtopos2(pairedcomps(sidx(j),3),:), EEG_chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp'); % full set, EEG
        title(['Comp ' num2str(pairedcomps(sidx(j),3)) ', p = ' num2str(pvalue(nulldistr(:), pairedcomps(sidx(j),5)), '%.3f')])
        topo2 = subplot(3, 10, meglocs2(i));
        topoplot(megtopos2(pairedcomps(sidx(j),3),:), MEG_chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp'); % full set, MEG
        axis([axis * 1.2]);
        clim1 = max(abs(caxis(topo1)));
        clim2 = max(abs(caxis(topo2)));
        clim = max(clim1, clim2);
        caxis(topo1, [-1 1]*clim);
        caxis(topo2, [-1 1]*clim);

        if i == 5 || j == size(pairedcomps,1)
            p = p+1;
            suptitle(['Split-half matched components for ' sublist{subno} ', sorted by midfrontalness: ' num2str(p*5 - 4) ' - ' num2str(p * 5)]);
            colormap bluewhitered;
            saveas(gcf, [dirs.plots 'r2r_' sublist{subno} '_' num2str(p) '_splithalves.png']);
            close;
        end
    end
end