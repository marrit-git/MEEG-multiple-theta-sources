%% Perform assorted supporting analyses, statistical tests, and ad hoc inspections.
% Not all of these analyses are reported in the associated manuscript; some
% are sanity checks, some are uninteresting.
%
% Analysis code for Simon task MEEG dataset.
%
% Author: Marrit Zuure
% 2019

%% Set paths
dirs = setpaths();

%% Set data import preliminaries
[sublist, ~, ~] = getICs2remove();

%% Loop over subjects
for subno = 1:length(sublist)
    disp(['Processing subject ' num2str(subno) ' of ' num2str(length(sublist)) '...']);
    
    %% Load GED data and previous analysis results
    GED_filename = [dirs.results sublist{subno} '_GED.mat'];
    load(GED_filename);
    
    ana_filename = [dirs.results sublist{subno} '_ana.mat'];
    load(ana_filename);
    
     %% Do some pre-computations for plots outside the subject loop
    % For average TF decomposition weighted by relative eigenvalue:
    % calculate relative component eigenvalues
    GED.rel_evals = GED.evals ./ sum(GED.evals(midf.comps2use));
    % Scale mean relative eigenvalue to 1
    GED.rel_evals = GED.rel_evals ./ mean(GED.rel_evals);
    
    %% Add data to structure array for pooling across subjects
    
     data(subno) =  struct('EEG', EEG, 'MEG', MEG, 'MEEG', MEEG, 'GED', GED, ...
        'tfd', tfd, 'thetacorr', thetacorr, 'compsynch', compsynch, 'trialtype', trialtype, ...
        'midf', midf, 'allrtidx', allrtidx, 'taskmod', taskmod, 'conflictmod', conflictmod, ...
        'sensor', sensor, 'elecsynch', elecsynch, 'gc', gc, 'noisecorr', noisecorr);
    
end

% Load PCA results (not subject-specific)
pca_filename = [dirs.results 'pca.mat'];
load(pca_filename);

%% Determine pairwise synchrony between components vs. electrodes

compsynch = [];
elecsynch = [];
for subno = 1:9
    compsynch = [compsynch, mean(mean(data(subno).compsynch.synch,2))];
    elecsynch = [elecsynch, mean(data(subno).elecsynch.synch)];
end
disp(['compsynch: ', num2str(compsynch), ' mean: ' num2str(mean(compsynch)), ' SD: ' num2str(std(compsynch))]);
disp(['elecsynch: ', num2str(elecsynch), ' mean: ' num2str(mean(elecsynch)), ' SD: ' num2str(std(elecsynch))]);

% Two-sided T-test
[h, p, ~, stats] = ttest2(compsynch, elecsynch, 'tail', 'both');
disp(['h = ', num2str(h), ', p = ' num2str(p)]);
disp(stats);

%% Compute component theta power time series pairwise correlations, both within and across trials

corr_within = [];
corr_across = [];
corr_across_cC = [];
corr_across_cI = [];
for subno = 1:9
    triu_idx = logical(triu(ones(data(subno).midf.num_comps),1));
    corr_within = [corr_within, mean(data(subno).thetacorr.within_trials(triu_idx))];
    corr_across = [corr_across, mean(data(subno).thetacorr.across_trials(triu_idx))];
    trialamp_cC = data(subno).thetacorr.trial_power(:,data(subno).trialtype==1);
    trialamp_cI = data(subno).thetacorr.trial_power(:,data(subno).trialtype==2);
    corr_cC = corr(trialamp_cC');
    corr_cI = corr(trialamp_cI');
    corr_across_cC = [corr_across_cC, mean(corr_cC(triu_idx))];
    corr_across_cI = [corr_across_cI, mean(corr_cI(triu_idx))];
end

disp(['corr_within: ', num2str(corr_within), ' mean: ' num2str(mean(corr_within)), ' SD: ' num2str(std(corr_within))]);
disp(['corr_across: ', num2str(corr_across), ' mean: ' num2str(mean(corr_across)), ' SD: ' num2str(std(corr_across))]);
disp(['corr_across_cC: ', num2str(corr_across_cC), ' mean: ' num2str(mean(corr_across_cC)), ' SD: ' num2str(std(corr_across_cC))]);
disp(['corr_across_cI: ', num2str(corr_across_cI), ' mean: ' num2str(mean(corr_across_cI)), ' SD: ' num2str(std(corr_across_cI))]);

% Check across-trial correlations being (significantly or not) *lower* than 0
[h,p, ~, stats] = ttest(corr_across, 0, 'tail', 'left'); % note that these are subject means, not means of pairs.
disp(['significantly negative across-trial correlations: h ' num2str(h), ' p ' num2str(p)]);
disp(stats);

%% Theta power deflections in TF decompositions

frex = dsearchn(data(1).tfd.frex(:), [5 7]'); % frequency window
widx = dsearchn(data(1).tfd.times2save(:), [0 800]'); % time window

thres = .5; % dB
comps = 0;
threscomps = 0;
for subno = 1:9
   comps = comps + data(subno).midf.num_comps;
   threscomps = threscomps + sum(max(max(squeeze(data(subno).tfd.tf(:, frex(1):frex(2), widx(1):widx(2), 1)), [], 3), [], 2) > thres);
end
disp([num2str(threscomps) '/' num2str(comps) ' comp peaks above ' num2str(thres) 'dB']);

%% Theta power component task and conflict modulation

taskmod = [];
confmod = [];
thetapwr_task = [];
thetapwr_baseline = [];
thetapwr_inc = [];
thetapwr_cong = [];
for subno = 1:9
    taskmod = [taskmod; data(subno).taskmod.score];
    confmod = [confmod; data(subno).conflictmod.score];
    
    widx = dsearchn(data(subno).tfd.times2save', data(subno).taskmod.window');
    bidx = dsearchn(data(subno).tfd.times2save', data(subno).tfd.baseline');

    thetaidx = dsearchn(data(subno).tfd.frex(:), data(subno).tfd.theta);
    
    thetapwr_task = [thetapwr_task; mean(data(subno).tfd.tf(:,thetaidx,widx(1):widx(2),1),3)];
    thetapwr_baseline = [thetapwr_baseline; mean(data(subno).tfd.tf(:,thetaidx,bidx(1):bidx(2),1),3)];
    
    thetapwr_inc = [thetapwr_inc; mean(data(subno).tfd.tf(:,thetaidx,widx(1):widx(2),3),3)];
    thetapwr_cong = [thetapwr_cong; mean(data(subno).tfd.tf(:,thetaidx,widx(1):widx(2),2),3)];
end

disp(['Component mean theta power task ' num2str(mean(thetapwr_task)) ', SD ' num2str(std(thetapwr_task))]);
disp(['Component mean theta power baseline ' num2str(mean(thetapwr_baseline)) ', SD ' num2str(std(thetapwr_baseline))]);
disp(['Component mean theta power congruent ' num2str(mean(thetapwr_cong)) ', SD ' num2str(std(thetapwr_cong))]);
disp(['Component mean theta power incongruent ' num2str(mean(thetapwr_inc)) ', SD ' num2str(std(thetapwr_inc))]);

% Test theta power on task against theta power on baseline
[h,p,~,stats] = ttest(thetapwr_task, thetapwr_baseline, 'tail', 'right');
disp(['Task modulation h: ' num2str(h) ', p: ' num2str(p) ', df = ' num2str(stats.df) ', t = ' num2str(stats.tstat) ]);

% Test theta power on incongruent against theta power on congruent trials
[h,p,~,stats] = ttest(thetapwr_inc, thetapwr_cong, 'tail', 'right');
disp(['Conflict modulation h: ' num2str(h) ', p: ' num2str(p) ', df = ' num2str(stats.df) ', t = ' num2str(stats.tstat) ]);

%% Theta power task and conflict modulation at FCz

fcz_taskmod = [];
fcz_confmod = [];

fcz_task = [];
fcz_baseline = [];
fcz_inc = [];
fcz_cong = [];
for subno = 1:9
    fczidx = strcmpi('fcz', {data(subno).MEEG.chanlocs.labels});
    fcz_taskmod = [fcz_taskmod, data(subno).sensor.task(fczidx)];
    fcz_confmod = [fcz_confmod, data(subno).sensor.conflict(fczidx)];
    
    widx = dsearchn(data(subno).tfd.times2save(:), data(subno).sensor.window');
    bidx = dsearchn(data(subno).tfd.times2save(:), data(subno).tfd.baseline');
    thetaidx = dsearchn(data(subno).tfd.frex(:), data(subno).tfd.theta);
    
    fcz_task = [fcz_task; squeeze(mean(data(subno).sensor.fcz(:,thetaidx,widx(1):widx(2),1),3))];
    fcz_baseline = [fcz_baseline; squeeze(mean(data(subno).sensor.fcz(:,thetaidx,bidx(1):bidx(2),1),3))];
    fcz_cong = [fcz_cong; squeeze(mean(data(subno).sensor.fcz(:,thetaidx,widx(1):widx(2),2),3))];
    fcz_inc = [fcz_inc; squeeze(mean(data(subno).sensor.fcz(:,thetaidx,widx(1):widx(2),3),3))];
end

disp(['FCz task mean ' num2str(mean(fcz_task)) ', SD ' num2str(std(fcz_task))]);
disp(['FCz baseline mean ' num2str(mean(fcz_baseline)) ', SD ' num2str(std(fcz_baseline))]);
disp(['FCz congruent mean ' num2str(mean(fcz_cong)) ', SD ' num2str(std(fcz_cong))]);
disp(['FCz incongruent mean ' num2str(mean(fcz_inc)) ', SD ' num2str(std(fcz_inc))]);

% Test FCz theta power on task against theta power on baseline
[h,p,~,stats] = ttest(fcz_task, fcz_baseline, 'tail', 'right');
disp(['FCz task modulation h: ' num2str(h) ', p: ' num2str(p) ', df = ' num2str(stats.df) ', t = ' num2str(stats.tstat) ]);

% Test FCz theta power on incongruent against theta power on congruent trials
[h,p,~,stats] = ttest(fcz_inc, fcz_cong, 'tail', 'right');
disp(['FCz conflict modulation h: ' num2str(h) ', p: ' num2str(p) ', df = ' num2str(stats.df) ', t = ' num2str(stats.tstat) ]);

%% GC component task and conflict modulation

dmtaskmod = [];
dmconfmod = [];
rmtaskmod = [];
rmconfmod = [];
for subno = 1:9
    dmtaskmod = [dmtaskmod; data(subno).gc.driver_task];
    dmconfmod = [dmconfmod; data(subno).gc.driver_conflict];
    rmtaskmod = [rmtaskmod; data(subno).gc.receiver_task];
    rmconfmod = [rmconfmod; data(subno).gc.receiver_conflict];
end

disp(['Component DM task mod mean ' num2str(mean(dmtaskmod)) ', SD ' num2str(std(dmtaskmod))]);
disp(['Component DM conflict mod mean ' num2str(mean(dmconfmod)) ', SD ' num2str(std(dmconfmod))]);
disp(['Component RM task mod mean ' num2str(mean(rmtaskmod)) ', SD ' num2str(std(rmtaskmod))]);
disp(['Component RM conflict mod mean ' num2str(mean(rmconfmod)) ', SD ' num2str(std(rmconfmod))]);

%% T-test GC task and conflict modulation

DM_task = [];
RM_task = [];
DM_baseline = [];
RM_baseline = [];
DM_inc = [];
DM_cong = [];
RM_inc = [];
RM_cong = [];

for subno = 1:9
    widx = dsearchn(data(subno).gc.winmid(:), data(subno).gc.wlim');
    bidx = dsearchn(data(subno).gc.winmid(:), data(subno).gc.baseline');
    
    DM_task = [DM_task; mean(data(subno).gc.F_driver(:,widx(1):widx(2)),2)];
    RM_task = [RM_task; mean(data(subno).gc.F_receiver(:,widx(1):widx(2)),2)];
    DM_baseline = [DM_baseline; mean(data(subno).gc.F_driver(:,bidx(1):bidx(2)),2)];
    RM_baseline = [RM_baseline; mean(data(subno).gc.F_receiver(:,bidx(1):bidx(2)),2)];
    
    DM_inc = [DM_inc; mean(data(subno).gc.F_driver_cI(:,widx(1):widx(2)),2)];
    RM_inc = [RM_inc; mean(data(subno).gc.F_receiver_cI(:,widx(1):widx(2)),2)];
    DM_cong = [DM_cong; mean(data(subno).gc.F_driver_cC(:,widx(1):widx(2)),2)];
    RM_cong = [RM_cong; mean(data(subno).gc.F_receiver_cC(:,widx(1):widx(2)),2)];
end

disp(['Component mean DM task ' num2str(mean(DM_task)) ', SD ' num2str(std(DM_task))]);
disp(['Component mean DM baseline ' num2str(mean(DM_baseline)) ', SD ' num2str(std(DM_baseline))]);
disp(['Component mean DM congruent ' num2str(mean(DM_cong)) ', SD ' num2str(std(DM_cong))]);
disp(['Component mean DM incongruent ' num2str(mean(DM_inc)) ', SD ' num2str(std(DM_inc))]);
 
disp(['Component mean RM task ' num2str(mean(RM_task)) ', SD ' num2str(std(RM_task))]);
disp(['Component mean RM baseline ' num2str(mean(RM_baseline)) ', SD ' num2str(std(RM_baseline))]);
disp(['Component mean RM congruent ' num2str(mean(RM_cong)) ', SD ' num2str(std(RM_cong))]);
disp(['Component mean RM incongruent ' num2str(mean(RM_inc)) ', SD ' num2str(std(RM_inc))]);

% Test GC on task against GC on baseline
[h,p,~,stats] = ttest(DM_task, DM_baseline, 'tail', 'both');
disp(['DM task modulation h: ' num2str(h) ', p: ' num2str(p) ', df = ' num2str(stats.df) ', t = ' num2str(stats.tstat) ]);
[h,p,~,stats] = ttest(RM_task, RM_baseline, 'tail', 'both');
disp(['RM task modulation h: ' num2str(h) ', p: ' num2str(p) ', df = ' num2str(stats.df) ', t = ' num2str(stats.tstat) ]);

% Test GC on incongruent against GC on congruent trials
[h,p,~,stats] = ttest(DM_inc, DM_cong, 'tail', 'both');
disp(['DM conflict modulation h: ' num2str(h) ', p: ' num2str(p) ', df = ' num2str(stats.df) ', t = ' num2str(stats.tstat) ]);
[h,p,~,stats] = ttest(RM_inc, RM_cong, 'tail', 'both');
disp(['RM conflict modulation h: ' num2str(h) ', p: ' num2str(p) ', df = ' num2str(stats.df) ', t = ' num2str(stats.tstat) ]);

%% T-test GC conflict modulation; timepoint by timepoint

dm_cC = [];
dm_cI = [];
for subno = 1:9
    dm_cC = [dm_cC; data(subno).gc.F_driver_cC];
    dm_cI = [dm_cI; data(subno).gc.F_driver_cI];
end

clear h p;
for i = 1:size(dm_cC,2)
    [h(i), p(i)] = ttest2(dm_cC(:,i), dm_cI(:,i), 'tail', 'left');
end

% Multiple comparisons correction, cluster based -> actually, no need; no
% single time point is significant

% Plot:
figure;
subplot(2,1,1);
plot(data(subno).gc.winmid, mean(dm_cC, 1));
hold on;
plot(data(subno).gc.winmid, mean(dm_cI, 1));
xlim([-500 1000]);
title('GC cC vs. cI')
subplot(2,1,2);
plot(data(subno).gc.winmid, p);
xlim([-500 1000]);
hold on;
ylim([0 1]);
plot([xlim], [0.05 0.05], 'k'); % alpha level without multiple comparisons correction
title('P-value');

%% T-test GC task modulation vs. baseline; timepoint by timepoint

dm_baseline_cC = [];
dm_task_cC = [];
dm_baseline_cI = [];
dm_task_cI = [];
% testing just for cI right now
bidx = dsearchn(data(1).gc.winmid(:), data(1).gc.baseline');
widx = dsearchn(data(1).gc.winmid(:), data(1).gc.window');
for subno = 1:9
    dm_baseline_cC = [dm_baseline_cC; data(subno).gc.F_driver_cC(:, bidx(1):bidx(2))];
    dm_baseline_cI = [dm_baseline_cI; data(subno).gc.F_driver_cI(:, bidx(1):bidx(2))];
    dm_task_cC = [dm_task_cC; data(subno).gc.F_driver_cC];
    dm_task_cI = [dm_task_cI; data(subno).gc.F_driver_cI];
end

clear h p;
for i = 1:size(dm_task_cC,2)
    [h_cC(i), p_cC(i)] = ttest2(dm_task_cC(:,i), mean(dm_baseline_cC,2), 'tail', 'both');
    [h_cI(i), p_cI(i)] = ttest2(dm_task_cI(:,i), mean(dm_baseline_cI,2), 'tail', 'both'); % note: size(dm_baseline_cC) is 52 x 9
end

% Plot:
figure;
subplot(2,1,1);
plot(data(subno).gc.winmid, mean(dm_task_cC, 1), 'b', 'DisplayName', 'cC');
hold on;
plot(data(subno).gc.winmid, mean(dm_task_cI, 1), 'r', 'DisplayName', 'cI');
xlim([-500 1000]);
title('GC task vs baseline')
legend show;
subplot(2,1,2);
plot(data(subno).gc.winmid, p_cC, 'b');
hold on;
plot(data(subno).gc.winmid, p_cI, 'r');
xlim([-500 1000]);
ylim([0 1]);
xax = xlim;
plot([xlim], [0.05 0.05], 'k'); % alpha level without multiple comparisons correction
plot([xlim], [0.0188 0.0188], 'Y'); % alpha level without multiple comparisons correction
title('P-value');

% Again, not significant at any time point in window of interest

%% Explained variance PCA and MEG PCs, factors

% Compute percentages explained variance
compts_expvar = pcan.compts_evals ./ sum(pcan.compts_evals) * 100;
receiver_expvar = pcan.receiver_evals ./ sum(pcan.receiver_evals) * 100;
driver_expvar = pcan.driver_evals ./ sum(pcan.driver_evals) * 100;
MEG_expvar = pcan.MEG_evals ./ sum(pcan.MEG_evals) * 100;

disp('compts:');
disp(compts_expvar(1:10));

disp('driver:');
disp(driver_expvar(1:10));

disp('receiver:');
disp(receiver_expvar(1:10));

disp('MEG:');
disp(MEG_expvar(1:10));

%% Granger causality order in ms

order = data(subno).gc.order * 1000 / (data(subno).MEEG.srate / data(subno).gc.ds_factor); % in ms
disp(order);

%% Granger causality significance

sig_driver_cC = [];
sig_driver_cI = [];
sig_receiver_cC = [];
sig_receiver_cI = [];

% Gather data
for subno = 1:9
    sig_driver_cC = [sig_driver_cC; data(subno).gc.sig_driver_cC];
    sig_driver_cI = [sig_driver_cI; data(subno).gc.sig_driver_cI];
    sig_receiver_cC = [sig_receiver_cC; data(subno).gc.sig_receiver_cC];
    sig_receiver_cI = [sig_receiver_cI; data(subno).gc.sig_receiver_cI];
end

% Check how many time points are significant for each of the components
disp(['Driver cC: ' num2str(sum(sig_driver_cC,2)') ' out of ' num2str(size(sig_driver_cC,2)) ' time points significant']);
disp(['Driver cI: ' num2str(sum(sig_driver_cI,2)') ' out of ' num2str(size(sig_driver_cI,2)) ' time points significant']);
disp(['Receiver cC: ' num2str(sum(sig_receiver_cC,2)') ' out of ' num2str(size(sig_receiver_cC,2)) ' time points significant']);
disp(['Receiver cI: ' num2str(sum(sig_receiver_cI,2)') ' out of ' num2str(size(sig_receiver_cI,2)) ' time points significant']);

%% Find peak times of theta power PCs

subno = 1;
plot(pcan.compts_evecs(:,1));
[peaks, idx] = max(pcan.compts_evecs(:,1:2));
window = data(subno).tfd.times2save;
widx = dsearchn(data(subno).tfd.times2save(:), pcan.window');
window = window(widx(1):widx(2));
comptspeaks = window(idx);
disp(comptspeaks);

%% Check noisy self-correlations: at which amount of noise is there no significant difference between samples?

% Gather data
for subno = 1:9
    noise = data(subno).noisecorr.noise;
    % Check: how much noise needs to be added before ha, hw becomes 0  in all
    % cases ( = there is no significant difference in correlations between
    % the original and between the noisy copies)?
    idxw = find(1 - data(subno).noisecorr.hw, 1);
    idxa = find(1 - data(subno).noisecorr.ha, 1);
    
    % If not statistically indistinguishable anywhere, set to max
    if isempty(idxw)
        idxw = length(noise);
    end
    if isempty(idxa)
        idxa = length(noise);
    end
    
    noise_same_within(subno) = noise(idxw);
    noise_same_across(subno) = noise(idxa);
    wstats(subno) = data(subno).noisecorr.wstats(idxw);
    astats(subno) = data(subno).noisecorr.astats(idxa);
end

disp(['Within-trial: noisy eigenvector self-copies need to have at least ' num2str(min(noise_same_within)*100) '% noise added to become as dissimilar as the original components.'])
disp(['Cross-trial: noisy eigenvector self-copies need to have at least ' num2str(min(noise_same_across)*100) '% noise added to become as dissimilar as the original components.'])
disp('(steps of 10%, noise is white noise with SD_noise = noiselevel * SD_eigenvector)');
disp('Test statistics within trials:');
disp([wstats.tstat]);
disp('df:');
disp([wstats.df]);
disp('Test statistics across trials:');
disp([astats.tstat]);
disp('df:');
disp([astats.df]);

%% Compute theta increase in dB at FCz for S08 (mostly copied from sensors.m)

% Sensor-level analysis task and baseline windows
sensor.prestim = [-500 -100]; % in ms
sensor.poststim = [100 800]; % in ms

% Independent components to remove from EEG/MEG: (previously determined)
EEGICs2remove = {
    'S08'     [ 1 3 4 6 14  ];
};
MEGICs2remove = {
     'S08'     [ 2 13 ];
};
sublist = EEGICs2remove(:,1);
subno = 1;

disp(['Processing subject ' num2str(subno) ' of ' num2str(length(sublist)) '...']);

% Load MEEG data
[~, ~, MEEG.data, ~, ~] = loadMEEG(sublist, subno, EEGICs2remove, MEGICs2remove, dirs.data, 'data');
    
% TF decompose FCz at all frequencies
S08_tf = tfdecomp(MEEG.data, MEEG, 2, 20, 40, 'means', tfd.times2save, tfd.baseline, trialtype);

% Extract max theta (5-7 Hz) dB increase at FCz
fczidx = strcmpi('fcz', {MEEG.chanlocs.labels});
frex = dsearchn(data(1).tfd.frex(:), [5 7]'); % frequency window
widx = dsearchn(data(1).tfd.times2save(:), [0 800]'); % time window
[peak, peakidx] = max(max(squeeze(S08_tf(fczidx,frex(1):frex(2), widx(1):widx(2), 1))));

disp(['Rejected subject has theta peak of ' num2str(peak) ' dB at electrode FCz.']);

%% Get correlations between component driving mass, receiving mass time courses

dm = [];
rm = [];
for subno = 1:9
    dm = [dm; data(subno).gc.F_driver];
    rm = [rm; data(subno).gc.F_receiver];
end
 [dmrmcorr, dmrmp] = corrcoef(dm', rm');
 disp(['Correlations between component driving and receiving mass: ' num2str(dmrmcorr(1,2)) ' (p: ' num2str(dmrmp(1,2)) ')']);
