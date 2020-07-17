function [ EEG, MEG, MEEG, allrtidx, trialtype ] = loadMEEG( sublist, subno, EEGICs2remove, MEGICs2remove, datadir, fields )
%LOADMEEG Loads EEG/MEG data, removes ICs, Z-scores and merges data.
%   LOADMEEG( sublist, subno, EEGICs2remove, MEGICs2remove, datadir ) 
%   loads EEG data (from a preprocessed .mat) and MEG data (from
%   condition-specific nut .mat files) from the folder with the subject name
%   (sublist{subno}) in datadir. Error trials are removed, EEGICs2remove and 
%   MEGICs2remove are removed, baseline is removed.
%   The data are de-meaned and Z-scored to convert them to the same scale.
%
%   FUNCTION INPUTS
%       sublist         Cell list of subject names
%       subno           Integer index of current subject
%       EEGICs2remove   Dictionary; subject names as keys, pre-computed 
%                       independent component indices to remove as values
%       MEGICs2remove   Dictionary; subject names as keys, pre-computed 
%                       independent component indices to remove as values
%       datadir         Path of directory to access to load subject data
%       fields          'all' to return full structs; 'data' to return only
%                       EEG.data, MEG.data, MEEG.data. Default is 'all'.
% 
%   FUNCTION OUTPUTS
%       EEG         EEGlab-like struct containing EEG data
%       MEG         EEGlab-like struct containing MEG data
%       MEEG        EEGlab-like struct containing EEG and MEG data
%       allrtidx    Array of trial reaction time indices
%       trialtype   Array of trial types (1 = congruent (cC), 2 = incongruent (cI))
%
%   Note: this function depends on EEGlab being in the MATLAB path.

%% Pre-check input validity
if ~strcmpi(fields, 'all') && ~strcmpi(fields, 'data')
    error(['Call to loadMEEG failed: ''' fields ''' invalid value for input parameter ''fields'' (try ALL or DATA)']);
end

%% load EEG data

clear global
load([ datadir sublist{subno} '/' sublist{subno} '_EEGpreproc.mat' ])
o=EEG; clear EEG; EEG=o; clear o; clear ALLEEG % weird, but unglobifies EEG

% remove components
EEG = eeg_checkset(EEG);
EEG = pop_subcomp(EEG,EEGICs2remove{subno,2});

% remove rejected trials from markers
markers(rejected_trials,:) = [];

% find error trials and remove from markers and from EEG data
errortrials = markers(:,4)==0;
markers(errortrials,:) = [];
EEG = pop_select(EEG,'notrial',find(errortrials));
EEG.icaact=[];
EEG.data = double(EEG.data);

% relabel marker values for convenience
markers(:,5) = (markers(:,1)<3);     % 1 = congruent
markers(:,5) = (markers(:,1)>2) + 1; % 2 = incongruent

trialtype = markers(:,5);

disp('Loaded EEG data.');

%% hard-coded EEG data corrections
%if subject is S06, swap electrodes fz and pz, fcz and cz
if strcmp(sublist{subno}, 'S06')
    temp = EEG.chanlocs(17);
    EEG.chanlocs(17) = EEG.chanlocs(19);
    EEG.chanlocs(19) = temp;
    
    temp = EEG.chanlocs(18);
    EEG.chanlocs(18) = EEG.chanlocs(31);
    EEG.chanlocs(31) = temp;
end

% if subject is S01, drop electrode f4
if strcmp(sublist{subno}, 'S01')
    EEG.data = [EEG.data(1:3,:,:); EEG.data(5:end,:,:)];
    EEG.chanlocs = [EEG.chanlocs(1:3), EEG.chanlocs(5:end)];
    EEG.nbchan = size(EEG.data,1);
end

%% load in MEG data

%load MEG ICA information
megica = load([ datadir 'MEGica/' sublist{subno} '_MEGica.mat' ]);

condition_list = {'cC';'cI'}; % means markers(:,1)<3 vs. markers(:,1)>2
for condi=1:length(condition_list)
    
    % find data from this condition and all runs
    filenames = dir([ datadir sublist{subno} '/nut_*_' condition_list{condi} '*.mat' ]);
    
    % now load data from all runs for this condition
    for filei=1:length(filenames)
        
        % load new nuts
        new_nuts = load([ datadir sublist{subno} '/' filenames(filei).name ]);
        
        % combine all runs into one big nuts (and rts)
        if filei==1 && condi==1
            allnuts = new_nuts.nuts;
            allrtidx  = new_nuts.rts;
        else
            allnuts.meg.data = cat(3,allnuts.meg.data,new_nuts.nuts.meg.data);
            allrtidx           = cat(2,allrtidx,new_nuts.rts);
        end
    end
end

% remove EEG-rejected trials
allnuts.meg.data(:,:,rejected_trials) = [];
allrtidx(:,rejected_trials) = [];

% remove error trials
allnuts.meg.data(:,:,errortrials) = [];
allrtidx(:,errortrials)             = [];
allrtidx(1,isnan(allrtidx(1,:))) = round(nanmean(allrtidx(1,:)));

% convert to MEG (eeglab format)
ddims = size(allnuts.meg.data); % time X chan X trial
MEG.data = zeros(ddims([2 1 3]));

for chani=1:ddims(2)
    MEG.data(chani,:,:) = 1e15 * double(squeeze(allnuts.meg.data(:,chani,:)));
end

% create MEG structure
MEG.times  = EEG.times;
MEG.trials = ddims(3);
MEG.pnts   = ddims(1);
MEG.nbchan = size(MEG.data,1);
MEG.srate  = allnuts.meg.srate;
MEG.xmin   = MEG.times(1)/1000;
MEG.xmax   = MEG.times(end)/1000;

if ~isfield(MEG, 'setname')
    MEG.setname = [];
end
if ~isfield(MEG, 'icawinv') % if no ICA data have been loaded
    MEG.icawinv = [];
    MEG.icaweights = [];
    MEG.icasphere = [];
    MEG.icaact = [];
else
    MEG.icawinv     = megica.EEG.icawinv;
    MEG.icasphere   = megica.EEG.icasphere;
    MEG.icaweights  = megica.EEG.icaweights;
    MEG.icachansind = megica.EEG.icachansind;
    
    % remove ICs from MEG
    c2k = setdiff(1:size(MEG.icaweights,1),MEGICs2remove{subno,2});
    first = MEG.icawinv(:,c2k)*MEG.icaweights(c2k,:)*MEG.icasphere;
    for triali=1:MEG.trials
        MEG.data(:,:,triali) = first * MEG.data(:,:,triali);
    end
end

% get channel names and locations
for chani=1:MEG.nbchan
    MEG.chanlocs(chani).labels = allnuts.meg.sensor_labels{chani};
    MEG.chanlocs(chani).X = allnuts.meg.sensorCoord(chani,1,2);
    MEG.chanlocs(chani).Y = allnuts.meg.sensorCoord(chani,2,2);
    MEG.chanlocs(chani).Z = allnuts.meg.sensorCoord(chani,3,2);
end
MEG = pop_chanedit(MEG,'convert','cart2all');

% baseline correct
MEG = pop_rmbase(MEG,[-200 0]);

clear tmp new_nuts allnuts

disp('Loaded MEG data.');

%% z-norm and combine MEG and EEG data into one big dataset

EEG_znormed = zeros(size(EEG.data));
MEG_znormed = zeros(size(MEG.data));

% z-norm the EEG
for triali=1:EEG.trials
    EEG_znormed(:,:,triali) = zscore(EEG.data(:,:,triali));
end
% z-norm the MEG
for triali=1:MEG.trials
    MEG_znormed(:,:,triali) = zscore(MEG.data(:,:,triali));
end

% collect into MEEG.data
MEEG.data = double( cat(1,EEG_znormed,MEG_znormed) );
MEEG.nbchan = size(MEEG.data,1);

% clear large variables
clear EEG_znormed
clear MEG_znormed

%% copy fields containing useful information to MEEG

assert(EEG.srate == MEG.srate); % If not, one or both should be resampled

MEEG.trials = EEG.trials;
MEEG.times = EEG.times;
MEEG.srate = EEG.srate;
MEEG.pnts = EEG.pnts;

MEEG.chanlocs = [EEG.chanlocs, MEG.chanlocs];
% In practice these are never plotted together, but it may be useful to be
% able to access them without having to load the EEG and MEG variables

%% Check whether to return full structs, or just data fields
if strcmpi(fields, 'data')
    EEG = EEG.data;
    MEG = MEG.data;
    MEEG = MEEG.data;
end
