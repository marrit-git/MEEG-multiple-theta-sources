%% Perform statistics on behavioral (response time) data
%
% Analysis code for Simon task MEEG dataset. 
%
% Author: Marrit Zuure

%% Calculate the Gratton effect from behavioral data
% Gratton effect: iI faster than cI

close all; clear;

%% Set paths
dirs = setpaths();

%% Set data import preliminaries
[sublist, ~, ~] = getICs2remove();

%% Get mean RT per subject per condition
% Also extract errors (to compute accuracy)
condmeans = zeros(length(sublist), 4);
condacc = zeros(length(sublist), 4);
for subno = 1:length(sublist)
    load(sprintf('%s/%s_behresults.mat', dirs.behavior, sublist{subno}));
    for i = 1:length(condition_labels)
        condmeans(subno, i) = mean(markers(markers(:,1)==i & markers(:,4)==1,2));
        condacc(subno,i) = mean(markers(markers(:,1)==i,4));
    end
end

% Make a bar plot for quick inspection
bar(mean(condmeans));
xticklabels(condition_labels);
% Figure is consistent with both Simon and Gratton effect

% Perform two-way ANOVA on accuracy
% column factor = current conflict
% row factor = previous conflict
nobs = 9; % one observation per subject
accmtx = [condacc(:,1), condacc(:,2);
          condacc(:,3), condacc(:,4)];
[p1, tbl1, stats1] = anova2(accmtx, nobs);

% Perform two-way ANOVA on RT
rtmtx = [condmeans(:,1), condmeans(:,2);
         condmeans(:,3), condmeans(:,4)];
[p2, tbl2, stats2 ] = anova2(rtmtx, nobs);

