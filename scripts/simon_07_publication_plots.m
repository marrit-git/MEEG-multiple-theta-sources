%% Create and save publication-quality plots.
% Analysis code for Simon task MEEG dataset. Makes use of Nathan Childress' 
% bluewhitered colormap, version 1, and Jonas' distributionPlot, 1.15; both 
% available from the Mathworks File Exchange.
% Author: Marrit Zuure
% March 2019

close all; clear;

%% Set paths
dirs = setpaths();

%% Set data import preliminaries
[sublist, EEGICs2remove, MEGICs2remove] = getICs2remove();

%% Set MATLAB plotting defaults

lw = 1.5;   % line width
msz = 8;    % marker size
fsz = 12;   % font size
alw = 0.75; % axes line width
tlen = [0.01 0.01]; % tick length

width = 8.5;  % width in cm 
height = 6; % height in cm
% Note: J Neurosci requests 8.5 cm, 11.6 cm, or 17.6 cm (1, 1.5, 2 columns) as max width

%% Propagate MATLAB plotting defaults
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultAxesFontSize',fsz);
set(0, 'defaultAxesLineWidth', alw);
set(0, 'defaultAxesTickLength', tlen);

% Set default renderer to vector format outputs. Note: doesn't support
% alpha channels.
set(0, 'defaultFigureRenderer', 'painters')

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','centimeters');
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

close; % close ghost figure

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

% Now start the actual plotting!
%% 1) Task, behavior, and sensor-level analyses
fig = 1;
figdir = [dirs.plots 'Figure ' num2str(fig)];
if ~exist(figdir, 'dir')
    mkdir(figdir);
end
%% 1a) Task design
% Done entirely in Illustrator
%% 1b) Simon and Gratton effects
subfig = 'b';

% Get mean RT and variance per subject per condition
condmeans = zeros(length(sublist), 4);
condvar = zeros(length(sublist), 4);
n = zeros(length(sublist), 4);
for subno = 1:length(sublist)
    load(sprintf('%sbehavior/%s_behresults.mat', dirs.data, sublist{subno}));
    for i = 1:length(condition_labels)
        condvar(subno,i) = var(markers(markers(:,1)==i & markers(:,4)==1,2));
        condmeans(subno, i) = mean(markers(markers(:,1)==i & markers(:,4)==1,2));
        n(subno, i) = numel(markers(markers(:,1)==i));
    end
end

% variance to SD
condsd = sqrt(mean(condvar));
% SD to SEM
condsem = condsd ./ sqrt(mean(n));

% error bars (SEM)
errhigh = condsem;
errlow = condsem;

ax = axgrid(4, 4, .5, .5, .5, .5, 1, 1, 'centimeters');

% Make a bar plot
axes(ax(1));
bar(mean(condmeans));
hold on;
errorbar(1:4, mean(condmeans), errlow, errhigh, 'LineWidth', 1.5, 'LineStyle', 'none', 'color', [0 0 0]);
xticklabels(condition_labels);
ylim([400 500]);
xlabel('Trial type');
ylabel('Reaction time (ms)');
set(gca, 'box', 'off');

set(gcf, 'PaperPosition', [0 0 50 50]);
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;
%% 1c) TF decomposition at FCz
subfig = 'c';
fcz_pooled_all = [];
fcz_pooled_conflict = [];
for subno = 1:length(sublist)
    fcz_pooled_all = cat(3, fcz_pooled_all, squeeze(data(subno).sensor.fcz(1,:,:,1)));
    fcz_pooled_conflict = cat(3, fcz_pooled_conflict, squeeze(data(subno).sensor.fcz(1,:,:,3)) - ...
         squeeze(data(subno).sensor.fcz(1,:,:,2)));
end

clim1 = 4;
clim2 = 1;

% Frequency window to outline
lo = pcan.lo;
hi = pcan.hi;

ax = axgrid(5, 14, .5, .5, 1, .5, 1, 2, 'centimeters');

axes(ax(1));
title('all trials');
contourf(data(subno).tfd.times2save, data(subno).tfd.frex, mean(fcz_pooled_all,3), 40, 'linecolor', 'none');
hold on;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('all trials');
xlim([-200 1000]);
caxis([-1 1] * clim1);
colorbar;
colormap bluewhitered;
% add TF window
plot([0 800], [lo lo], 'k', 'LineWidth', 1);
plot([0 800], [hi hi], 'k', 'LineWidth', 1);
plot([0 0], [lo hi], 'k', 'LineWidth', 1);
plot([800 800], [lo hi], 'k', 'LineWidth', 1);

axes(ax(2));
contourf(data(subno).tfd.times2save, data(subno).tfd.frex, mean(fcz_pooled_conflict,3), 40, 'linecolor', 'none');
hold on;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('conflict modulation');
xlim([-200 1000]);
caxis([-1 1] * clim2);
colorbar;
colormap bluewhitered;
% add TF window
plot([0 800], [lo lo], 'k', 'LineWidth', 1);
plot([0 800], [hi hi], 'k', 'LineWidth', 1);
plot([0 0], [lo hi], 'k', 'LineWidth', 1);
plot([800 800], [lo hi], 'k', 'LineWidth', 1);

saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;
%% 1d) EEG and MEG task and conflict modulation topoplots
subfig = 'd';

% Sensor locations aren't identical acrosss subjects. Feed data through 
% meanffms to get a joint set of chanlocs and a matching mean 
% filter forward model back. Do EEG and MEG separately to avoid mixing
% EEG and MEG channels in the output.
for subno = 1:length(sublist)
   task_ffms_EEG{subno} = data(subno).sensor.task(1:data(subno).EEG.nbchan);
   task_ffms_MEG{subno} = data(subno).sensor.task(data(subno).EEG.nbchan+1:end);
   conflict_ffms_EEG{subno} = data(subno).sensor.conflict(1:data(subno).EEG.nbchan);
   conflict_ffms_MEG{subno} = data(subno).sensor.conflict(data(subno).EEG.nbchan+1:end);
   chanlocs_EEG{subno} = data(subno).EEG.chanlocs;
   chanlocs_MEG{subno} = data(subno).MEG.chanlocs;
end
[sensor_task_EEG, sensor_chanlocs_EEG_1] = meanffms(task_ffms_EEG, chanlocs_EEG, 'labels', length(sublist)-2, 2);
[sensor_task_MEG, sensor_chanlocs_MEG_1] = meanffms(task_ffms_MEG, chanlocs_MEG, 'labels', length(sublist)-2, 2);
[sensor_conflict_EEG, sensor_chanlocs_EEG_2] = meanffms(conflict_ffms_EEG, chanlocs_EEG, 'labels', length(sublist)-2, 2);
[sensor_conflict_MEG, sensor_chanlocs_MEG_2] = meanffms(conflict_ffms_MEG, chanlocs_MEG, 'labels', length(sublist)-2, 2);

clim1 = 2;
clim2 = 0.5;
clim3 = 2;
clim4 = 0.5;

figure;
subplot(2,2,1);
topoplot(sensor_task_EEG, sensor_chanlocs_EEG_1, 'electrodes', 'ptslabels', 'style', 'map', 'shading', 'interp');
colorbar;
caxis([-1 1]*clim1);
subplot(2,2,2);
topoplot(sensor_conflict_EEG, sensor_chanlocs_EEG_2, 'electrodes', 'ptslabels', 'style', 'map', 'shading', 'interp');
colorbar;
caxis([-1 1]*clim2);
subplot(2,2,3);
topoplot(sensor_task_MEG, sensor_chanlocs_MEG_1, 'style', 'map', 'shading', 'interp');
colorbar;
caxis([-1 1]*clim3);
axis([axis * 1.2]);
subplot(2,2,4);
topoplot(sensor_conflict_MEG, sensor_chanlocs_MEG_2, 'style', 'map', 'shading', 'interp');
colorbar;
caxis([-1 1]*clim4);
axis([axis * 1.2]);
colormap bluewhitered;

saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 2) Explanation of GED and data that goes into it
fig = 2;
figdir = [dirs.plots 'Figure ' num2str(fig)];
if ~exist(figdir, 'dir')
    mkdir(figdir);
end
%% 2a1) Filtered EEG/MEG sensor time series
subfig = 'a1';

% Pick a subject for sensor time series
subno = 3;

% Load MEEG data
disp('Loading data...');
[~, ~, MEEG.data] = loadMEEG(sublist, subno, EEGICs2remove, MEGICs2remove, dirs.data, 'data');
disp('Data loaded.');

% Pick a trial to show
trial = 19;

% Pick a few channels to show
EEGchannels = [11, 28, 31]; % under 56
MEGchannels = [30, 110, 221] + 56; % 56-328
channels = [EEGchannels, MEGchannels];

% Get broadband sensor data
bbdata = squeeze(MEEG.data(channels, :, trial));

% Bandpass filter the data using the original data analysis settings
disp('Filtering data...');
kernel = data(subno).GED.kernel;
fft_len = size(bbdata,2)+length(kernel)-1; % number of time points/frequencies for fft to return
trim_len = (length(kernel)-1)/2; % number of time points to trim from start and end of result
fdata = 2*real( ifft( bsxfun(@times,fft(bbdata,fft_len,2),fft(kernel, fft_len)) ,[],2) );
fdata = reshape(fdata(:,trim_len+1:end-trim_len,:), size(bbdata,1), size(bbdata,2));

% Trim broadband and filtered data to 0-800 ms window
widx = dsearchn(MEEG.times(:), [0 800]');
bbdata = bbdata(:, widx(1):widx(2));
fdata = fdata(:, widx(1):widx(2));

% Add offsets
distance = 2;
bbscale = 0.4;
for i = 1:size(bbdata,1)
    bbdata(i,:) = bbdata(i,:) * bbscale + i * distance;
    fdata(i,:) = fdata(i,:) + i * distance;
end

% Plot filtered trials
figure;
hold on;
for i = 1:3
    plot(fdata(i,:), 'b');
end
for i = 4:6
    plot(fdata(i,:), 'r');
end
set(gca, 'YDir','reverse'); % Plot EEG above, MEG below
axis off;
text(-120, 2*distance, 'EEG', 'FontSize', 20);
text(-120, 5*distance, 'MEG', 'FontSize', 20);
title('\theta-filtered')

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;
%% 2a2) Broadband EEG/MEG sensor time series
subfig = 'a2';

% Plot broadband trials
figure;
hold on;
for i = 1:3
    plot(bbdata(i,:), 'b');
end
for i = 4:6
    plot(bbdata(i,:), 'r');
end
set(gca, 'YDir','reverse'); % Plot EEG above, MEG below
axis off;
title('broadband');
text(-120, 2*distance, 'EEG', 'FontSize', 20);
text(-120, 5*distance, 'MEG', 'FontSize', 20);

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;
%% 2b) R, S, W, and Lambda for equation (and part of 2a)
subfig = 'b';

% Pick a subject for scree plot and topoplot
subno = 3;

clim_S = 1;
clim_W = .5;
clim_R = 1;
clim_Lambda = max(data(subno).GED.orig_evals) * .2;

figure;
h(1) = subplot(1,5,1);
imagesc(data(subno).GED.covS); % S
colormap(h(1), bluewhitered);
caxis([-1 1] * clim_S);
axis square;
xticks([]);
yticks([]);
xlabel('S')
h(2) = subplot(1,5,2);
imagesc(data(subno).GED.orig_evecs); % W
colormap(h(2), bluewhitered);
caxis([-1 1] * clim_W);
xticks([]);
yticks([]);
axis square;
xlabel('W')
h(3) = subplot(1,5,3);
imagesc(data(subno).GED.covR); % R
colormap(h(3), bluewhitered);
caxis([-1 1] * clim_R);
xticks([]);
yticks([]);
axis square;
xlabel('R')
h(4) = subplot(1,5,4);
imagesc(data(subno).GED.orig_evecs); % W
colormap(h(4), bluewhitered);
caxis([-1 1] * clim_W);
xticks([]);
yticks([]);
axis square;
xlabel('W')
h(5) = subplot(1,5,5);
imagesc(diag(data(subno).GED.orig_evals)); % Lambda
colormap(h(5), bluewhitered);
caxis([0 1] * clim_Lambda);
xticks([]);
yticks([]);
axis square;
xlabel('\Lambda')

% Set size
width = 17;
height = 10;
set(gcf, 'PaperPosition', [0 0 width height]);

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;
%% 2c) Illustration of component time series creation
%% 2c1) Raw single-trial data
subfig = 'c1';

comp = 1;
trial = 225;
sensors = 30:32;
subno = 3;

if ~isfield(MEEG, 'data')
    % Load MEEG data
    disp('Loading data...');
    [~, ~, MEEG.data] = loadMEEG(sublist, subno, EEGICs2remove, MEGICs2remove, dirs.data, 'data');
    disp('Data loaded.');
end

% get and plot sensors on trial
trialdata = MEEG.data(sensors,:,trial);
% Add offsets
distance = 2;
scale = 0.5;
for i = 1:size(trialdata,1)
    trialdata(i,:) = trialdata(i,:) * scale + i * distance;
end
figure;
plot(MEEG.times, trialdata);
xlim([-300 1000]);

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;
%% 2c2) Start of component eigenvector
subfig = 'c2';
% Get and plot component weights 1-3 (start of eigenvector)
evec = data(subno).GED.evecs(:, data(subno).midf.comps2use(comp));
evec = evec(sensors);
% Actually, just extract the color values and plot in Illustrator
image(evec*250);
colormap bluewhitered;

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;
%% 2c3) Component single-trial time course
subfig = 'c3';

% Get and plot component trial
comptrial = data(subno).midf.compts(comp,:,trial);
figure;
plot(MEEG.times, comptrial);
xlim([-300 1000]);

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 2d1) Scree plot of eigenvalues and significance threshold
subfig = 'd1';

evals_to_plot = sort([data(subno).GED.evals(:); data(subno).GED.repeat_evals(:)], 'descend');
figure;
h1 = plot(evals_to_plot, '.', 'MarkerSize', 20);
% set(h1, 'markerfacecolor', get(h1, 'color'));
hold on;
if all(data(subno).GED.evals_pval) == 0
    plot(data(subno).GED.perm_evals, '--k');
end
% legend('eigenvalues', 'significance threshold');
xlabel('Components');
ylabel('S to R power ratio (\lambda)');
ylim([-.15 max(evals_to_plot)+0.1]);
xlim([-5 328]); % hard-coding number of sensors
set(gca, 'Box', 'off');
xticks([]);
yticks([]);

% Set size
width = 17;
height = 10;
set(gcf, 'PaperPosition', [0 0 width height]);

% Plot a circle around the second eigenvalue
hold on;
x = 2;
y = evals_to_plot(2);
r = 6; % on x-axis
yscale = diff(ylim) / diff(xlim) * (width/height);
th = linspace(0, 2*pi, round(4 * pi * r)); % Define angles
x = r * cos(th) + x;
y = yscale * r * sin(th) + y;
plot(x, y, 'r', 'LineWidth', 2, 'HandleVisibility', 'off');
text(200, 1.1, 'p = 0.05', 'Color', 'r');

% Set size again; fixes crop in export somehow
set(gcf, 'PaperPosition', [0 0 width height]);

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;
%% 2d2) Component time course and EEG and MEG topography for component from scree plot
subfig = 'd2';
figure;
subplot(2,2,[1 2]);
plot(mean(data(subno).GED.compts(2,widx(1):widx(2),:),3)); % mean over trials
axis off;
topo1 = subplot(2,2,3);
topoplot(data(subno).midf.ffm_EEG(:,2), data(subno).EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
topo2 = subplot(2,2,4);
topoplot(data(subno).midf.ffm_MEG(:,2), data(subno).MEG.chanlocs,  'electrodes', 'off', 'style', 'map', 'shading', 'interp');
axis([axis * 1.2]);
colormap(bluewhitered);

% Use same color axis for both plots
clim1 = max(abs(caxis(topo1)));
clim2 = max(abs(caxis(topo2)));
clim = max(clim1, clim2);
caxis(topo1, [-1 1]*clim);
caxis(topo2, [-1 1]*clim);

% Set size
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 3) Component retention steps + uniqueness scatterplots + top 2 comps per subject
fig = 3;
figdir = [dirs.plots 'Figure ' num2str(fig)];
if ~exist(figdir, 'dir')
    mkdir(figdir);
end
%% 3a) EEG template; EEG + MEG topo 1; EEG + MEG topo 2 (from same subject)
subfig = 'a';

% Choose subject
subno = 3;
comps = [1, 6];

figure;
subplot(3,2,1);
topoplot(data(subno).midf.template, data(subno).EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
subplot(3,2,3);
topoplot(data(subno).midf.ffm_EEG(:,comps(1)), data(subno).EEG.chanlocs, 'electrodes', 'off', 'style', 'map',  'shading', 'interp');
subplot(3,2,4);
topoplot(data(subno).midf.ffm_MEG(:,comps(1)), data(subno).MEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
axis([axis * 1.2]);
subplot(3,2,5);
topoplot(data(subno).midf.ffm_EEG(:,comps(2)), data(subno).EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
subplot(3,2,6);
topoplot(data(subno).midf.ffm_MEG(:,comps(2)), data(subno).MEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
axis([axis * 1.2]);
colormap(bluewhitered);

% % Set size
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;
%% 3b) 9 times Scree plot zoom, significance bounding box, color-coded markers
subfig = 'b';

figure;
% set up figure axes
ax(1) = axes('Position', [.04375, .6813, .275, .275]);
ax(2) = axes('Position', [.3625, .6813, .275, .275]);
ax(3) = axes('Position', [.6813, .6813, .275, .275]);
ax(4) = axes('Position', [.04375, .3625, .275, .275]);
ax(5) = axes('Position', [.3625, .3625, .275, .275]);
ax(6) = axes('Position', [.6813, .3625, .275, .275]);
ax(7) = axes('Position', [.04375, .04375, .275, .275]);
ax(8) = axes('Position', [.3625, .04375, .275, .275]);
ax(9) = axes('Position', [.6813, .04375, .275, .275]);

for subno = 1:length(sublist)
    axes(ax(subno));
    evals_to_plot_orig = sort([data(subno).GED.evals(:); data(subno).GED.repeat_evals(:)], 'descend');
    evals_idx = 1:length(evals_to_plot_orig);
   
    % Remove evals_to_plot that are in repeat_evals or in
    % midf.comps2use (don't want double markers)
    evals_to_plot = evals_to_plot_orig;
    evals_to_plot(data(subno).GED.repeat_evals_idx) = [];
    evals_idx(data(subno).GED.repeat_evals_idx) = [];
    evals_to_plot(data(subno).midf.comps2use) = [];
    evals_idx(data(subno).midf.comps2use) = [];
    
    plot(evals_idx, evals_to_plot, '.', 'DisplayName', 'components', 'MarkerSize', 20,  'color', [.6 .6 .6]);
    hold on;
    
    % Fix midfrontal marker offsets by dropped repeat eigenvalues (were
    % previously overwriting blue markers)
    cum_repeat = cumsum(data(subno).GED.repeat_evals_idx); % ugh @ variable name
    midf_idx = data(subno).midf.comps2use + cum_repeat(data(subno).midf.comps2use);
    
    plot(find(data(subno).GED.repeat_evals_idx), data(subno).GED.repeat_evals, '.r', 'DisplayName', '\Delta\lambda < 1%', 'MarkerSize', 20); % Dropped: in red
    plot(midf_idx, evals_to_plot_orig(midf_idx), '.g', 'DisplayName', 'R^2 > 0.5', 'MarkerSize', 20); % Midfrontal: in green

    plot(data(subno).GED.perm_evals, '--k', 'DisplayName', 'p = 0.05');
    
    plotbottom = data(subno).GED.perm_evals(find(evals_to_plot_orig < data(subno).GED.perm_evals', 1)); % plot perm_evals dashed line instead of bounding box
    
%     Temporarily disabled to prevent Matlab from stretching plot:
%     title(['Subject ' num2str(subno)]);
%     xlabel('Component');
%     ylabel('Power ratio (\lambda)')
%     legend show;
    ylim([plotbottom - .1, max(evals_to_plot_orig)+0.1]);
    xlim([-1 sum(evals_to_plot_orig > data(subno).GED.perm_evals')+5]);
    box off;
    
end

% Set size
width = 21;
height = 15;
set(gcf, 'PaperPosition', [0 0 width height]);

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 3c) EEG/MEG topographies per subject + subject average + PCA
subfig = 'c';

% This plot contains many vertices and takes a lot of memory to generate
% and save. Splitting it into multiple figures and recombining in
% Illustrator.

figure;
ax = axgrid(15, 7, .2, .2, .2, .2, length(sublist)+2, 4, 'centimeters');
plotno = 1:4:length(sublist)*4;
for subno = 1:3
    % EEG topographies comp 1
    axes(ax(plotno(subno)));
    topo1 = ax(plotno(subno));
    comps2use = data(subno).midf.comps2use;
    topoplot(data(subno).midf.ffm_EEG(:,comps2use(1)), data(subno).EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
    % MEG topographies comp 1
    axes(ax(plotno(subno)+1));
    topo2 = ax(plotno(subno)+1);
    comps2use = data(subno).midf.comps2use;
    topoplot(data(subno).midf.ffm_MEG(:,comps2use(1)), data(subno).MEG.chanlocs, 'style', 'map', 'shading', 'interp');
    axis([axis * 1.2]);
    
    % EEG topographies comp 2
    axes(ax(plotno(subno)+2));
    topo1 = ax(plotno(subno)+2);
    comps2use = data(subno).midf.comps2use;
    topoplot(data(subno).midf.ffm_EEG(:,comps2use(2)), data(subno).EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
    % MEG topographies comp 2
    axes(ax(plotno(subno)+3));
    topo2 = ax(plotno(subno)+3);
    comps2use = data(subno).midf.comps2use;
    topoplot(data(subno).midf.ffm_MEG(:,comps2use(2)), data(subno).MEG.chanlocs, 'style', 'map', 'shading', 'interp');
    axis([axis * 1.2]);
end
colormap bluewhitered;
saveas(gcf, sprintf('%sFigure %i/%i%s-1.pdf', dirs.plots, fig, fig, subfig));
close;

figure;
ax = axgrid(15, 7, .2, .2, .2, .2, length(sublist)+2, 4, 'centimeters');
for subno = 4:6
    % EEG topographies comp 1
    axes(ax(plotno(subno)));
    topo1 = ax(plotno(subno));
    comps2use = data(subno).midf.comps2use;
    topoplot(data(subno).midf.ffm_EEG(:,comps2use(1)), data(subno).EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
    % MEG topographies comp 1
    axes(ax(plotno(subno)+1));
    topo2 = ax(plotno(subno)+1);
    comps2use = data(subno).midf.comps2use;
    topoplot(data(subno).midf.ffm_MEG(:,comps2use(1)), data(subno).MEG.chanlocs, 'style', 'map', 'shading', 'interp');
    axis([axis * 1.2]);
    
    % EEG topographies comp 2
    axes(ax(plotno(subno)+2));
    topo1 = ax(plotno(subno)+2);
    comps2use = data(subno).midf.comps2use;
    topoplot(data(subno).midf.ffm_EEG(:,comps2use(2)), data(subno).EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
    % MEG topographies comp 2
    axes(ax(plotno(subno)+3));
    topo2 = ax(plotno(subno)+3);
    comps2use = data(subno).midf.comps2use;
    topoplot(data(subno).midf.ffm_MEG(:,comps2use(2)), data(subno).MEG.chanlocs, 'style', 'map', 'shading', 'interp');
    axis([axis * 1.2]);
end
colormap bluewhitered;
saveas(gcf, sprintf('%sFigure %i/%i%s-2.pdf', dirs.plots, fig, fig, subfig));
close;

figure;
ax = axgrid(15, 7, .2, .2, .2, .2, length(sublist)+2, 4, 'centimeters');
for subno = 7:length(sublist)
    % EEG topographies comp 1
    axes(ax(plotno(subno)));
    topo1 = ax(plotno(subno));
    comps2use = data(subno).midf.comps2use;
    topoplot(data(subno).midf.ffm_EEG(:,comps2use(1)), data(subno).EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
    % MEG topographies comp 1
    axes(ax(plotno(subno)+1));
    topo2 = ax(plotno(subno)+1);
    comps2use = data(subno).midf.comps2use;
    topoplot(data(subno).midf.ffm_MEG(:,comps2use(1)), data(subno).MEG.chanlocs, 'style', 'map', 'shading', 'interp');
    axis([axis * 1.2]);
    
    % EEG topographies comp 2
    axes(ax(plotno(subno)+2));
    topo1 = ax(plotno(subno)+2);
    comps2use = data(subno).midf.comps2use;
    topoplot(data(subno).midf.ffm_EEG(:,comps2use(2)), data(subno).EEG.chanlocs, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
    % MEG topographies comp 2
    axes(ax(plotno(subno)+3));
    topo2 = ax(plotno(subno)+3);
    comps2use = data(subno).midf.comps2use;
    topoplot(data(subno).midf.ffm_MEG(:,comps2use(2)), data(subno).MEG.chanlocs, 'style', 'map', 'shading', 'interp');
    axis([axis * 1.2]);
end
colormap bluewhitered;
saveas(gcf, sprintf('%sFigure %i/%i%s-3.pdf', dirs.plots, fig, fig, subfig));
close;

% Create subject average topoplots
figure;
ax = axgrid(15, 7, .2, .2, .2, .2, length(sublist)+2, 4, 'centimeters');
for subno = 1:length(sublist)
    EEG1{subno} = data(subno).midf.ffm_EEG(:,comps2use(1));
    MEG1{subno} = data(subno).midf.ffm_MEG(:,comps2use(1));
    EEG2{subno} = data(subno).midf.ffm_EEG(:,comps2use(2));
    MEG2{subno} = data(subno).midf.ffm_MEG(:,comps2use(2));
    EEGchan{subno} = data(subno).EEG.chanlocs;
    MEGchan{subno} = data(subno).MEG.chanlocs;
end
[EEG1_mean, newEEGchan] = meanffms(EEG1, EEGchan, 'labels');
[MEG1_mean, newMEGchan] = meanffms(MEG1, MEGchan, 'labels');
[EEG2_mean, ~] = meanffms(EEG2, EEGchan, 'labels');
[MEG2_mean, ~] = meanffms(MEG2, MEGchan, 'labels');

% Plot subject average topoplots
axes(ax(length(sublist)*4+1));
topo1 = ax(length(sublist)*4+1);
topoplot(EEG1_mean, newEEGchan, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
axes(ax(length(sublist)*4+2));
topo2 = ax(length(sublist)*4+2);
topoplot(MEG1_mean, newMEGchan, 'style', 'map', 'shading', 'interp');
axis([axis * 1.2]);

axes(ax(length(sublist)*4+3));
topo1 = ax(length(sublist)*4+3);
topoplot(EEG2_mean, newEEGchan, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
axes(ax(length(sublist)*4+4));
topo2 = ax(length(sublist)*4+4);
topoplot(MEG2_mean, newMEGchan, 'style', 'map', 'shading', 'interp');
axis([axis * 1.2]);
colormap bluewhitered;
saveas(gcf, sprintf('%sFigure %i/%i%s-4.pdf', dirs.plots, fig, fig, subfig));
close;

% Create topoplot PCs
subfig = 'd';

% Gather EEG and apply PCA
figure;
ax = axgrid(15, 7, .2, .2, .2, .2, length(sublist)+2, 4, 'centimeters');
EEG_ffms = {};
EEG_chanlocs = {};
for subno = 1:length(sublist)
    EEG_ffms{subno} = data(subno).midf.ffm_EEG;
    EEG_chanlocs{subno} = data(subno).EEG.chanlocs;
end
[~, newEEGchan, newEEGffms] = meanffms(EEG_ffms, EEG_chanlocs, 'labels');
% Clear rows with NaNs before computing covariance matrix
nonans = logical(1-any(isnan(newEEGffms')));
newEEGffms = newEEGffms(nonans,:);
newEEGchan = newEEGchan(nonans);
% Calculate covariance matrix
EEG_cov = newEEGffms * newEEGffms';
% PCA
[EEG_evecs, EEG_evals] = eig(EEG_cov);
% Sort by eigenvalue
[EEG_evals, sidx] = sort(diag(EEG_evals), 'descend');
EEG_evecs = EEG_evecs(:, sidx);

% Gather MEG and apply PCA
MEG_ffms = {};
MEG_chanlocs = {};
for subno = 1:length(sublist)
    MEG_ffms{subno} = data(subno).midf.ffm_MEG;
    MEG_chanlocs{subno} = data(subno).MEG.chanlocs;
end
[~, newMEGchan, newMEGffms] = meanffms(MEG_ffms, MEG_chanlocs, 'labels');
% Calculate covariance matrix
MEG_cov = newMEGffms * newMEGffms';
% PCA
[MEG_evecs, MEG_evals] = eig(MEG_cov);
% Sort by eigenvalue
[MEG_evals, sidx] = sort(diag(MEG_evals), 'descend');
MEG_evecs = MEG_evecs(:, sidx);

% Create MEG and EEG PC topoplots
EEGPC1 = EEG_evecs(:,1)' * EEG_cov;
EEGPC2 = EEG_evecs(:,2)' * EEG_cov;
MEGPC1 = MEG_evecs(:,1)' * MEG_cov;
MEGPC2 = MEG_evecs(:,2)' * MEG_cov;

% Plot topoplot PCs
% subplot(length(sublist)+2,4,length(sublist)*4+5);
axes(ax(length(sublist)*4+5));
topoplot(-1 * EEGPC1, newEEGchan, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
% subplot(length(sublist)+2,4,length(sublist)*4+6);
axes(ax(length(sublist)*4+6));
topoplot(EEGPC2, newEEGchan, 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
% subplot(length(sublist)+2,4,length(sublist)*4+7);
axes(ax(length(sublist)*4+7));
topoplot(MEGPC1, newMEGchan, 'style', 'map', 'shading', 'interp');
axis([axis * 1.2]);
% subplot(length(sublist)+2,4,length(sublist)*4+8);
axes(ax(length(sublist)*4+8));
topoplot(MEGPC2, newMEGchan, 'style', 'map', 'shading', 'interp');
axis([axis * 1.2]);
colormap bluewhitered;

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 4) Component TF decomposition averages
fig = 4;
figdir = [dirs.plots 'Figure ' num2str(fig)];
if ~exist(figdir, 'dir')
    mkdir(figdir);
end

subfig = 'a';
% Extract values and multiply component TF decompositions by relative component eigenvalues
tfd_pooled_all = [];
tfd_pooled_conflict = [];
for subno = 1:length(sublist)
    scale = data(subno).GED.rel_evals(1:data(subno).midf.num_comps);
    tfd_pooled_all = cat(1, tfd_pooled_all, squeeze(data(subno).tfd.tf(:,:,:,1)) .* scale);
    tfd_pooled_conflict = cat(1, tfd_pooled_conflict, (squeeze(data(subno).tfd.tf(:,:,:,3)) - squeeze(data(subno).tfd.tf(:,:,:,2))) .* scale);
end

% Take mean
mean_tfd_all = squeeze(mean(tfd_pooled_all, 1));
mean_tfd_conflict = squeeze(mean(tfd_pooled_conflict, 1));

clim1 = 17;
clim2 = 3;

ax = axgrid(5, 14, .5, .5, 1, .5, 1, 2, 'centimeters');

axes(ax(1));
title('all trials');
contourf(data(subno).tfd.times2save, data(subno).tfd.frex, mean_tfd_all, 40, 'linecolor', 'none');
hold on;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('all trials');
xlim([-200 1000]);
caxis([-1 1] * clim1);
colorbar;
colormap bluewhitered;

axes(ax(2));
contourf(data(subno).tfd.times2save, data(subno).tfd.frex, mean_tfd_conflict, 40, 'linecolor', 'none');
hold on;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('conflict modulation');
xlim([-200 1000]);
caxis([-1 1] * clim2);
colorbar;
colormap bluewhitered;

saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 5) Pairwise synchrony + correlation scatterplots
fig = 5;
figdir = [dirs.plots 'Figure ' num2str(fig)];
if ~exist(figdir, 'dir')
    mkdir(figdir);
end
subfig = 'abc';

% Gather data: for each subject, lowest noise level at which self-correlations
% are significantly different from component pairwise correlations
iterations = data(1).noisecorr.iterations;
noise = data(1).noisecorr.noise;

noise_diff_within = zeros(length(sublist));
noise_diff_across = zeros(length(sublist));
noise_distr_within = cell(length(sublist),1);
noise_distr_across = cell(length(sublist),1);
for subno = 1:length(sublist)
    % within trials
    firstdiff = find(1 - data(subno).noisecorr.hw, 1); % find first h = 0 for significant difference in correlations
    if isempty(firstdiff)
        firstdiff_w(subno) = length(data(subno).noisecorr.hw);
    else
        firstdiff_w(subno) = firstdiff;
    end
    
    % across trials
    firstdiff = find(1 - data(subno).noisecorr.ha, 1); % find first h = 0 for significant difference in correlations
    if isempty(firstdiff)
        firstdiff_a(subno) = length(data(subno).noisecorr.ha);
    else
        firstdiff_a(subno) = firstdiff;
    end
end
% extract minimum across subjects
firstdiff_w = min(firstdiff_w);
firstdiff_a = min(firstdiff_a);

for subno = 1:length(sublist)
    triu_idx = logical(triu(ones(iterations, iterations),1));
    noise_distr_within{subno} = squeeze(data(subno).noisecorr.within(:,firstdiff_w,triu_idx)); % all self-correlation values for all components
    noise_distr_across{subno} = squeeze(data(subno).noisecorr.across(:,firstdiff_a,triu_idx)); % all self-correlation values for all components
end
% Reshape each components x self-correlations matrix into a single vector
% of all self-correlations
for subno = 1:length(sublist)
    noise_distr_across{subno} = reshape(noise_distr_across{subno}, 1, []);
    noise_distr_within{subno} = reshape(noise_distr_within{subno}, 1, []);
end

% Collect data
pwsynch = [];
subno_pairs = []; % same pairs of components per subject for multiple plots, only need to construct this once
elecsynch = [];
for subno = 1:length(sublist)
    triu_idx = find(triu(data(subno).compsynch.mtx,1));
    pwsynch = [pwsynch; data(subno).compsynch.mtx(triu_idx)];
    subno_pairs = [subno_pairs; subno * ones(1,length(triu_idx))'];
    elecsynch = [elecsynch, mean(data(subno).elecsynch.synch)];
end

withincorr = [];
acrosscorr = [];
for subno = 1:length(sublist)
    triu_idx = find(triu(data(subno).thetacorr.within_trials,1));
    withincorr = [withincorr; data(subno).thetacorr.within_trials(triu_idx)];
    acrosscorr = [acrosscorr; data(subno).thetacorr.across_trials(triu_idx)];
    triu_idx = find(triu(squeeze(data(subno).noisecorr.within(1,1,:,:)), 1));
end

clear ax;
figure;
% set up figure axes
% margins:
ml = 1.5;
mr = .1;
mb = 1.5;
mt = .1;
rows = 1;
cols = 3;
ax = axgrid(5, 17.6, mb, mt, ml, mr, rows, cols, 'centimeters');

axes(ax(1));
colors = get(gca, 'ColorOrder');
scatter(subno_pairs, pwsynch, 'filled'); % adding alpha back in later in Illustrator
hold on;
scatter(1:length(sublist), elecsynch, 'r', 'filled'); % Pz and POz pairwise synchrony
xlabel('Subject');
ylabel('Pairwise synchrony');
xlim([0.1 9.9]); %hard-coding number of subjects
xticks(1:9);
ylim([-0.08*0.6 0.6]);
legend('Component pairs', 'Electrodes Pz-POz')
set(gca, 'box', 'off');
set(gca,'FontSize',9);

axes(ax(2));
distributionPlot(noise_distr_within', 'color', 'b', 'showMM', 6);
hold on;
scatter(subno_pairs, withincorr, [], colors(1,:), 'filled');
xlabel('Subject');
ylabel('Within-trial theta amplitude correlations');
xlim([0.1 9.9]); %hard-coding number of subjects
xticks(1:9);
ylim([-0.08 0.45]); % max is 0.21
% legend('Component pairs', 'Noisy self-correlations')
set(gca, 'box', 'off');
set(gca,'FontSize',9);

axes(ax(3));
distributionPlot(noise_distr_across', 'color', 'b', 'showMM', 6);
hold on;
scatter(subno_pairs, acrosscorr, [], colors(1,:), 'filled');
xlabel('Subject');
ylabel('Across-trial theta amplitude correlations');
xlim([0.1 9.9]); %hard-coding number of subjects
xticks(1:9);
ylim([-0.08 0.45]); % max is 0.40
% legend('Component pairs', 'Noisy self-correlations')
set(gca, 'box', 'off');
set(gca,'FontSize',9);

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 6) Theta power, task and conflict modulation scatterplots
fig = 6;
figdir = [dirs.plots 'Figure ' num2str(fig)];
if ~exist(figdir, 'dir')
    mkdir(figdir);
end

%% 6a) Conflict modulation
subfig = 'a';

% Create grand average, dB-normalized power
for subno = 1:length(sublist)
   thetaidx = dsearchn(data(subno).tfd.frex(:), data(subno).tfd.theta);
   cCpwr(subno,:) = mean(data(subno).tfd.tf(:,thetaidx,:,2),1); % mean over components for each subject, cC trials
   cIpwr(subno,:) = mean(data(subno).tfd.tf(:,thetaidx,:,3),1); % mean over components for each subject, cI trials
end
grandavg_cC = mean(cCpwr);
grandavg_cI = mean(cIpwr);

% Normalize per condition
bidx = dsearchn(data(1).tfd.times2save', data(1).tfd.baseline');
grandavg_cC = grandavg_cC - mean(grandavg_cC(bidx(1):bidx(2)));
grandavg_cI = grandavg_cI - mean(grandavg_cI(bidx(1):bidx(2)));

thetafreq = dsearchn(data(subno).tfd.frex(:), 6);
window = data(subno).conflictmod.window;

ax = axgrid(3.5, 4.5, .5, .5, .5, .5, 1, 1, 'centimeters');

axes(ax(1));
hold on;
plot(data(subno).tfd.times2save, grandavg_cC); % congruent
plot(data(subno).tfd.times2save, grandavg_cI); % incongruent
axy = ylim;
xlabel('Time (ms)');
ylabel('\theta power (dB)')
plot([window(1) window(1)], [axy(1) axy(2)], '--k', 'HandleVisibility','off'); % vertical lines for window
plot([window(2) window(2)], [axy(1) axy(2)], '--k', 'HandleVisibility','off');
xlim([-600 1000]);
ylim([-.5 2]);
set(gca, 'box', 'off');

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 6b) Task modulation scores for components and FCz
subfig = 'b';
% Collect data
taskmod = [];
subno_comps = []; % same number of components for every subject, so only needs to be generated once
fcz_taskmod = [];
for subno = 1:length(sublist)
    taskmod = [taskmod; data(subno).taskmod.score];
    subno_comps = [subno_comps, subno * ones(1, length(data(subno).taskmod.score))];
    fczidx = strcmpi('fcz', {data(subno).MEEG.chanlocs.labels});
    fcz_taskmod = [fcz_taskmod, data(subno).sensor.task(fczidx)];
end

ax = axgrid(3.5, 4, .5, .5, .5, .5, 1, 1, 'centimeters');

axes(ax(1));
scatter(subno_comps, taskmod, 'filled', 'MarkerFaceAlpha', .3, 'MarkerEdgeAlpha', .3);
hold on;
scatter(1:length(sublist), fcz_taskmod, 'r', 'filled');
% legend('Components', 'FCz');
xlabel('Subject');
ylabel('\theta power task modulation');
xticks(1:9);
xlim([0.1 9.9]); %hard-coding number of subjects
set(gca, 'box', 'off');

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 6c) Conflict modulation scores for components and FCz
subfig = 'c';

% Collect data
conflictmod = [];
fcz_conflictmod = [];
for subno = 1:length(sublist)
    conflictmod = [conflictmod; data(subno).conflictmod.score];
    fczidx = strcmpi('fcz', {data(subno).MEEG.chanlocs.labels});
    fcz_conflictmod = [fcz_conflictmod, data(subno).sensor.conflict(fczidx)];
end

ax = axgrid(3.5, 4, .5, .5, .5, .5, 1, 1, 'centimeters');

axes(ax(1));
scatter(subno_comps, conflictmod, 'filled', 'MarkerFaceAlpha', .3, 'MarkerEdgeAlpha', .3);
hold on;
scatter(1:length(sublist), fcz_conflictmod, 'r', 'filled');
% legend('Components', 'FCz'); % Commented out to prevent plot stretching
xlabel('Subject');
ylabel('\mu_inc - \mu_cong (dB)');
xticks(1:9);
xlim([0.1 9.9]); %hard-coding number of subjects
set(gca, 'box', 'off');

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 6d) Component time series heat maps
subfig = 'd';

freqidx = dsearchn(data(subno).tfd.frex(:), [pcan.lo pcan.hi]'); % mean over theta frequency range
compts = cell2mat( arrayfun(@(c) squeeze(mean(c(:).tfd.tf(:,freqidx(1):freqidx(2),:,1), 2))', data, 'Uniform', 0));
widx = dsearchn(data(subno).tfd.times2save(:), pcan.window');

numcomps = zeros(length(sublist),1);
for subno = 1:length(sublist)
    numcomps(subno) = sum(subno_comps == subno);
end
cumcomps = cumsum(numcomps);

ax = axgrid(12, 5.75 + 1, .5, .5, .5, .5, 3, 1, 'centimeters');

% Component time series
axes(ax(1));
imagesc(compts');
hold on;
for i = 1:length(sublist)-1
    plot([0 length(data(subno).tfd.times2save)], [cumcomps(i)+.5 cumcomps(i)+.5], 'k', 'LineWidth', .5);
end
xticks(dsearchn(data(subno).tfd.times2save(:), [-500:500:1000]'));
xticklabels(-500:500:1000);
xlim([min(xticks), max(xticks)]);
xlabel('Time (ms)');
caxis([-5 5]);
colormap(bluewhitered);
colorbar;

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 6e) Top PC line plots
subfig = 'e';

freqidx = dsearchn(data(subno).tfd.frex(:), [pcan.lo pcan.hi]'); % mean over 4-7 Hz
compts = cell2mat( arrayfun(@(c) squeeze(mean(c(:).tfd.tf(:,freqidx(1):freqidx(2),:,1), 2))', data, 'Uniform', 0));

widx = dsearchn(data(subno).tfd.times2save(:), pcan.window');

ax = axgrid(10, 4.75, .5, .5, .5, .5, 3, 1, 'centimeters');

axes(ax(1));
for i = 1:pcan.compts_PCs
    plot(data(subno).tfd.times2save(widx(1):widx(2)), pcan.compts_evecs(:,i), 'DisplayName', ['principal component ' num2str(i)]);
    hold on;
    xlim(pcan.window);
    xlabel('Time (ms)');
    ylabel('a.u.')
    set(gca, 'box', 'off');
end

axes(ax(2));
for i = 1:pcan.driver_PCs
    plot(data(subno).gc.winmid(widx_gc(1):widx_gc(2)), pcan.driver_evecs(:,i), 'DisplayName', ['principal component ' num2str(i)]);
    hold on;
    xlim(pcan.window);
    xlabel('Time (ms)');
    ylabel('a.u.')
    set(gca, 'box', 'off');
end

axes(ax(3));
for i = 1:pcan.receiver_PCs
    plot(data(subno).gc.winmid(widx_gc(1):widx_gc(2)), pcan.receiver_evecs(:,i), 'DisplayName', ['principal component ' num2str(i)]);
    hold on;
    xlim(pcan.window);
    xlabel('Time (ms)');
    ylabel('a.u.')
    set(gca, 'box', 'off');
end

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 7) Granger causality calculation; GC task and conflict modulation scatterplots
fig = 7;
figdir = [dirs.plots 'Figure ' num2str(fig)];
if ~exist(figdir, 'dir')
    mkdir(figdir);
end
%% 7a) Component time series to illustrate GC calculation
subfig = 'a';

subno = 3;
comps = 1:3;
trials = 105:107;
compts = data(subno).midf.compts(comps, :, trials);

ax = axgrid(2.5, 6, .5, .5, .2, .2, 1, length(comps), 'centimeters');

distance = .03;
scale = .3;
for j = 1:length(comps)
    axes(ax(j));
    hold on;
    for i = 1:length(trials)
        plot(MEEG.times, compts(j,:,i) * scale + (1+length(trials)-i)*distance);
    end
    ylim([0 0.2]);
    xlim([0 500]) % not the window used for GC calculations, but for illustrative purposes
    
    % Add marker for 200 ms window
    widx = dsearchn(MEEG.times(:), [150 350]');
    plot([MEEG.times(widx(1)) MEEG.times(widx(2))], [distance distance]);

    % Add marker for 50 ms time step
    widx = dsearchn(MEEG.times(:), [100 150]');
    plot([MEEG.times(widx(1)) MEEG.times(widx(2))], [distance distance]);
    
    box off;
end

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 7b1) GC DM conflict modulation line plot
subfig = 'b1';

% Create grand averages
for subno = 1:length(sublist)
    dm_cC(subno,:) = mean(data(subno).gc.F_driver_cC ./ data(subno).gc.driver_baseline, 1);
    dm_cI(subno,:) = mean(data(subno).gc.F_driver_cI ./ data(subno).gc.driver_baseline, 1);
end
grandavg_dm_cC = mean(dm_cC, 1);
grandavg_dm_cI = mean(dm_cI, 1);

window = data(subno).gc.window;

ax = axgrid(3.5, 4, .5, .5, .5, .5, 1, 1, 'centimeters');

axes(ax(1));
colors = get(gca,'ColorOrder');
plot(data(subno).gc.winmid, grandavg_dm_cC);
hold on;
plot(data(subno).gc.winmid, grandavg_dm_cI);
xlim([-600 1000]);
ylim([.85 1.15]);
axy = ylim;
xlabel('Time (ms)');
ylabel('GC driving mass');
plot([window(1) window(1)], [axy(1) axy(2)], '--k', 'HandleVisibility','off');
plot([window(2) window(2)], [axy(1) axy(2)], '--k', 'HandleVisibility','off');
xticks([0 800]);
set(gca, 'box', 'off');

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 7b2) GC RM conflict modulation line plot
subfig = 'b2';

% Create grand averages
for subno = 1:length(sublist)
    rm_cC(subno,:) = mean(data(subno).gc.F_receiver_cC ./ data(subno).gc.receiver_baseline, 1);
    rm_cI(subno,:) = mean(data(subno).gc.F_receiver_cI ./ data(subno).gc.receiver_baseline, 1);
end
grandavg_rm_cC = mean(rm_cC, 1);
grandavg_rm_cI = mean(rm_cI, 1);

window = data(subno).gc.window;

ax = axgrid(3.5, 4, .5, .5, .5, .5, 1, 1, 'centimeters');

axes(ax(1));
colors = get(gca,'ColorOrder');
plot(data(subno).gc.winmid, grandavg_rm_cC, 'color', colors(1,:));
hold on;
plot(data(subno).gc.winmid, grandavg_rm_cI, 'color', colors(2,:));
xlim([-600 1000]);
ylim([.85 1.15]);
axy = ylim;
xlabel('Time (ms)');
ylabel('GC receiving mass');
plot([window(1) window(1)], [axy(1) axy(2)], '--k', 'HandleVisibility','off');
plot([window(2) window(2)], [axy(1) axy(2)], '--k', 'HandleVisibility','off');
xticks([0 800]);
set(gca, 'box', 'off');

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;
%% 7c1) GC task modulation scatterplot
subfig = 'c1';

driver_taskmod_pooled = [];
receiver_taskmod_pooled = [];
for subno = 1:length(sublist)
    driver_taskmod_pooled = [driver_taskmod_pooled; data(subno).gc.driver_task];
    receiver_taskmod_pooled = [receiver_taskmod_pooled; data(subno).gc.receiver_task];
end

ax = axgrid(3.5, 4, .5, .5, .5, .5, 1, 1, 'centimeters');

axes(ax(1));
scatter(subno_comps-.2, driver_taskmod_pooled, 10, 'filled', 'MarkerFaceAlpha', .3, 'MarkerEdgeAlpha', .3);
hold on;
scatter(subno_comps+.2, receiver_taskmod_pooled, 10, 'filled', 'MarkerFaceAlpha', .3, 'MarkerEdgeAlpha', .3);
xlabel('Subject');
ylabel('\mu_{win} / \mu_{base}');
xlim([0.1 9.9]); %hard-coding number of subjects
xticks(1:9);
set(gca, 'box', 'off');
legend('driver', 'receiver');

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;
%% 7c2) GC conflict modulation scatterplot
subfig = 'c2';

ax = axgrid(3.5, 4, .5, .5, .5, .5, 1, 1, 'centimeters');

driver_conflictmod_pooled = [];
receiver_conflictmod_pooled = [];
for subno = 1:length(sublist)
    driver_conflictmod_pooled = [driver_conflictmod_pooled; data(subno).gc.driver_conflict];
    receiver_conflictmod_pooled = [receiver_conflictmod_pooled; data(subno).gc.receiver_conflict];
end

axes(ax(1));
scatter(subno_comps-.2, driver_conflictmod_pooled, 10, 'filled', 'MarkerFaceAlpha', .3, 'MarkerEdgeAlpha', .3);
hold on;
scatter(subno_comps+.2, receiver_conflictmod_pooled, 10, 'filled', 'MarkerFaceAlpha', .3, 'MarkerEdgeAlpha', .3);
xlabel('Subject');
ylabel('\mu_{inc} - \mu_{cong}');
xlim([0.1 9.9]); %hard-coding number of subjects
xticks(1:9);
set(gca, 'box', 'off');
legend('driver', 'receiver');

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 7d) Granger time series heat maps
subfig = 'd';

driverts = cell2mat( arrayfun(@(c) c(:).gc.F_driver', data, 'Uniform', 0));
receiverts = cell2mat( arrayfun(@(c) c(:).gc.F_receiver', data, 'Uniform', 0));
widx_gc = dsearchn(data(subno).gc.winmid(:), pcan.window');

numcomps = zeros(length(sublist),1);
for subno = 1:length(sublist)
    numcomps(subno) = sum(subno_comps == subno);
end
cumcomps = cumsum(numcomps);

ax = axgrid(12, 5.75 + 1, .5, .5, .5, .5, 3, 1, 'centimeters');

% Driver time series
axes(ax(2));
imagesc(driverts' ./ mean(driverts',2) - 1); % Normalized for display purposes, centering around 0
hold on;
for i = 1:length(sublist)-1
    plot([0 length(data(subno).tfd.times2save)], [cumcomps(i)+.5 cumcomps(i)+.5], 'k', 'LineWidth', .5);
end
xticks(dsearchn(data(subno).gc.winmid(:), [-500:500:1000]'));
xticklabels(-500:500:1000);
xlim([min(xticks), max(xticks)]);
caxis([-1.5 1.5]);
colormap(bluewhitered);
colorbar;

% Receiver time series
axes(ax(3));
imagesc(receiverts' ./ mean(receiverts',2) - 1); % Normalized for display purposes, centering around 0
hold on;
for i = 1:length(sublist)-1
    plot([0 length(data(subno).tfd.times2save)], [cumcomps(i)+.5 cumcomps(i)+.5], 'k', 'LineWidth', .5);
end
xticks(dsearchn(data(subno).gc.winmid(:), [-500:500:1000]'));
xticklabels(-500:500:1000);
xlim([min(xticks), max(xticks)]);
caxis([-1.5 1.5]);
colormap(bluewhitered);
colorbar;

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%% 7e) Top PC line plots
subfig = 'e';

freqidx = dsearchn(data(subno).tfd.frex(:), [pcan.lo pcan.hi]'); % mean over 4-7 Hz
compts = cell2mat( arrayfun(@(c) squeeze(mean(c(:).tfd.tf(:,freqidx(1):freqidx(2),:,1), 2))', data, 'Uniform', 0));

widx = dsearchn(data(subno).tfd.times2save(:), pcan.window');

ax = axgrid(10, 4.75, .5, .5, .5, .5, 3, 1, 'centimeters');

axes(ax(2));
for i = 1:pcan.driver_PCs
    plot(data(subno).gc.winmid(widx_gc(1):widx_gc(2)), pcan.driver_evecs(:,i), 'DisplayName', ['principal component ' num2str(i)]);
    hold on;
    xlim(pcan.window);
    xlabel('Time (ms)');
    ylabel('a.u.')
    set(gca, 'box', 'off');
end

axes(ax(3));
for i = 1:pcan.receiver_PCs
    plot(data(subno).gc.winmid(widx_gc(1):widx_gc(2)), pcan.receiver_evecs(:,i), 'DisplayName', ['principal component ' num2str(i)]);
    hold on;
    xlim(pcan.window);
    xlabel('Time (ms)');
    ylabel('a.u.')
    set(gca, 'box', 'off');
end

% Save
saveas(gcf, sprintf('%sFigure %i/%i%s.pdf', dirs.plots, fig, fig, subfig));
close;

%%
disp('Finished plotting.');
quit();