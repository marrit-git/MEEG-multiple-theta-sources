function dirs = setpaths()
% SETPATHS specifies data paths, output paths, EEGlab and Brainstorm installation paths, etc. 
% Edit the paths in this function before running the analysis code on a new machine.

% Input/output folder paths
dirs.home = '/home/mzuure/Data/Simon MEEG/'; % parent directory

dirs.scripts = [dirs.home 'scripts/']; % analysis scripts
dirs.data = [dirs.home 'data/']; % cleaned EEG and MEG data
dirs.behavior = [dirs.data 'behavior/']; % behavioral data
dirs.results = [dirs.home 'results/']; % to output/load intermediate analysis results
dirs.plots = [dirs.home 'figures/']; % to output figures

% Directories containing installations and helper packages (i.e.,
% everything to be added to path).
% NOTE: These are not included with the Simon MEEG manuscript code.
% Download them yourself and point these paths towards them.
pkgdirs.eeglab = '/home/mzuure/Documents/eeglab/'; % EEGlab, needed for many analyses
pkgdirs.mvgc = '/home/mzuure/MVGCtoolbox/'; % Multivariate Granger Causality toolbox, needed for GC analyses
pkgdirs.bluewhitered = [dirs.home 'packages/bluewhitered'];
pkgdirs.violinplot = [dirs.home '/packages/distributionPlot'];

% Create directories if they don't already exist
warning('off', 'MATLAB:MKDIR:DirectoryExists');
fn = fieldnames(dirs);
for i = 1:length(fn)
    status = mkdir(dirs.(fn{i}));
end

% Add package directories to MATLAB search path
fn = fieldnames(pkgdirs);
for i = 1:length(fn)
    addpath(genpath(pkgdirs.(fn{i})));
end
addpath(genpath(pkgdirs.eeglab));
addpath(genpath(pkgdirs.mvgc));

disp('Successfully set paths.');

end