function [ mean_ffm, new_chanlocs, new_ffms ] = meanffms( ffms, chanlocs, joinby, varargin )
%MEANFFMS( FFMS, CHANLOCS, JOINBY, MIN_OBSERVATIONS, FIRST2USE ) deals with taking 
% the mean over the filter forward model (topoplots), given that channel 
% locations may vary across subjects.
%   MEANFFMS extracts all unique coordinate sets from CHANLOCS, recombines
%   them into a new set of chanlocs, and returns it together with the mean
%   FFM for that set of chanlocs. To prevent channels with just a few
%   observations from dominating the mean in that location, a minimum
%   number of observations can be specified; any channel with fewer
%   observations than that will have all observations set to NaN. Note that
%   this function is designed to work with CHANLOCS containing XYZ coordinates.
%   Make sure that all passed FFMs are sign-flipped so that taking the mean
%   makes sense.
%   Other notes: passing joint EEG + MEG chanlocs might result in the
%   channels being unordered in the output, making plotting difficult. 
%   When joining by labels, note that the coordinates will be taken from
%   first2use (if available); this function does not calculate the mean
%   coordinates for a given label.
%   MEANFFMS requires EEGlab to be in the Matlab path.
%   MEANFFMS requires the extractfield() function from the Signal
%   Processing toolbox, which unfortunately can't be distributed along with
%   the rest of this code.
%   EEGlab version used: 14.1.2, released May 17th, 2018
%
%   FUNCTION INPUTS
%       ffms                (subjects x 1) cell array, containing (n x 
%                               channels) ffms over which to take the mean
%       chanlocs            (subjects x 1) cell array, containing EEGlab-
%                               style channel locations for each subject
%       joinby              What to use to determine channel identity: 
%                           'labels' or 'coordinates'.
%       min_observations    Optional: Integer value; minimum number of 
%                               observations per channel. Channels with 
%                               fewer are NaN'd. Defaults to 0.
%       first2use           Optional: Integer value; which subjects' 
%                               chanlocs to use first to match coordinates
%                               and labels. Useful if some subjects have 
%                               wrong electrode labels. Defaults to 1.
%
%   FUNCTION OUTPUTS
%       mean_ffm            Mean over input FFMs, matching new_chanlocs
%       new_chanlocs        Recombined set of all channel coordinates
%       new_ffms            Input FFMs, reordered to match new_chanlocs

%% Set default argument values
if isempty(varargin)
    min_observations = 0;
    first2use = 1;
elseif length(varargin) == 1
    min_observations = varargin{1};
    first2use = 1;
elseif length(varargin) > 1
    min_observations = varargin{1};
    first2use = varargin{2};
end

%% Extract unique coordinates and create new set of chanlocs

% Set reference chanlocs
ref_chanlocs = chanlocs{first2use};

% Set up initial coordinate and labels matrix
X = extractfield(ref_chanlocs, 'X');
Y = extractfield(ref_chanlocs, 'Y');
Z = extractfield(ref_chanlocs, 'Z');
labels = extractfield(ref_chanlocs, 'labels');

% Further populate coordinates and labels
for subno = 1:length(ffms)
    if subno == first2use
        continue; % skip; this subject has already been processed
    end
    
    X = [X, extractfield(chanlocs{subno}, 'X')];
    Y = [Y, extractfield(chanlocs{subno}, 'Y')];
    Z = [Z, extractfield(chanlocs{subno}, 'Z')];
    labels = [labels, extractfield(chanlocs{subno}, 'labels')];
end

% Join X, Y, Z into XYZ so unique combinations can be extracted
XYZ = [X; Y; Z];

if strcmpi(joinby, 'coordinates')
    % Extract unique combinations of coordinates; extract accompanying labels
    [unique_XYZ, idx, ~] = unique(XYZ', 'rows', 'stable');
    labels = labels(idx);
    
    % Create new chanlocs (struct array) containing unique coordinates
    for i = 1:size(unique_XYZ,1)
        new_chanlocs(i).labels = labels{i};
        new_chanlocs(i).X = unique_XYZ(i,1);
        new_chanlocs(i).Y = unique_XYZ(i,2);
        new_chanlocs(i).Z = unique_XYZ(i,3);
    end
elseif strcmpi(joinby, 'labels')
    % Extract unique labels; extract accompanying coordinates
    [unique_labels, idx, ~] = unique(labels, 'stable');
    XYZ = XYZ(:, idx);
    
    % Create new chanlocs (struct array) containing unique labels
    for i = 1:size(unique_labels,2)
        new_chanlocs(i).labels = unique_labels{i};
        new_chanlocs(i).X = XYZ(1,i);
        new_chanlocs(i).Y = XYZ(2,i);
        new_chanlocs(i).Z = XYZ(3,i);
    end
end

% Use EEGlab to re-populate all non-XYZ coordinates
new_chanlocs = pop_chanedit(new_chanlocs, 'convert', 'cart2all');

%% Rearrange all FFMs to match new_chanlocs, based on channel labels           
            
% Get total number of FFMs
total_ffms = 0;
for subno = 1:length(ffms)
    total_ffms = total_ffms + size(ffms{subno},2);
end

% Move all FFMs to the right place in new_ffms, which matches
% new_chanlocs
new_ffms = nan(length(new_chanlocs),total_ffms);
ffms_before = 0;
sum_observations = zeros(1,length(new_chanlocs));
for subno = 1:length(ffms)
    for i = 1:size(ffms{subno},1)
        move2chan = strcmpi(chanlocs{subno}(i).labels, {new_chanlocs.labels});
        new_ffms(move2chan, ffms_before+1:ffms_before+size(ffms{subno},2)) = ffms{subno}(i,:);
        sum_observations = sum_observations + move2chan;
    end
    ffms_before = ffms_before + size(ffms{subno},2);
end

%% If too few observations for a channel, set all observations to NaN

if exist('min_observations', 'var')
    new_ffms(sum_observations < min_observations, :) = nan;
end

%% Take mean FFM
mean_ffm = nanmean(new_ffms,2);

end

