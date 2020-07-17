%% Analyze differences in eigenvector modality between significant and 
% non-significant components.
% Also perform a control analysis to see if GED respects data unimodality.
% Analysis code for Simon task MEEG dataset, newly added during revision.
% Author: Marrit Zuure
% June 2019

close all; clear;

dirs = setpaths();

% Load EEGlab
eeglab;

[sublist, EEGICs2remove, MEGICs2remove] = getICs2remove();

%% Load data for all subjects

% Initialize variables
em = {};
md = {};
num_comps = zeros(length(sublist),1);

for subno = 1:length(sublist)
    disp(['Loading data for subject ' num2str(subno) '...']);
    load([dirs.results sublist{subno} '_GED.mat']); % Load file locally

    % Compute modality dominance value for GED eigenvectors (with
    % repeating eigenvalues removed)
    for i = 1:size(GED.evecs,2)
        e = rms(GED.evecs(1:EEG.nbchan,i));
        m = rms(GED.evecs(EEG.nbchan+1:end,i));
        em{subno}(i) = (e-m)/(e+m); % this value indicates whether EEG or MEG is dominant
        md{subno}(i) = ((e-m)/(e+m))^2; % this is the magnitude of the dominance ( = unimodality)
    end
    
    num_comps(subno) = GED.num_comps;
end

%% Compute overall MD stats (comparing significant vs.non-significant components)

md_sign = [];
md_nonsign = [];
for subno = 1:length(sublist)
    md_sign = [md_sign, md{subno}(1:num_comps(subno))];
    md_nonsign = [md_nonsign, md{subno}(num_comps(subno)+1:end)];
end

disp(mean(md_sign));
disp(std(md_sign));
disp(mean(md_nonsign));
disp(std(md_nonsign));

%T-test: significant comps md different from non-significant comps?
%(pooled)
[h_pooled,p_pooled,~,stats] = ttest2(md_sign, md_nonsign, 'tail', 'left');
disp('pooled:');
disp(['H = ' num2str(h_pooled)]);
disp(['p = ' num2str(p_pooled)]);
disp(stats);

%T-test: significant comps md different from non-significant comps?
%(per subject)
for subno = 1:length(sublist)
    nc_sig(subno) = num_comps(subno);
    nc_nonsig(subno) = length(md{subno});
    md_sign = md{subno}(1:num_comps(subno));
    md_nonsign = md{subno}(num_comps(subno)+1:end);
    [h,p,~,stats] = ttest2(md_sign, md_nonsign, 'tail', 'left');
    H(subno) = h;
    P(subno) = p;
    Stats{subno} = stats;
    disp([sublist{subno} ':']);
    disp(['H = ' num2str(h)]);
    disp(['p = ' num2str(p)]);
    disp(stats);
end
% Benjamini-Hochberg multiple comparisons correction:
all_p = [P, p_pooled];
label = {sublist{:}, 'pooled'};
[all_p, sidx] = sort(all_p, 'ascend');
for i = 1:length(all_p)
    disp([label{sidx(i)} '   p: ' num2str(all_p(i)) '   p*num_tests/rank: ' num2str(all_p(i)*length(all_p)/i) '   signif.: ' num2str((all_p(i)*length(all_p)/i) < 0.05)]);
end

% Report numbers of significant/non-significant comps per subject to
% justify pooled analysis
disp(mean(nc_sig))
disp(mean(nc_nonsig))

% Next: control analysis: does GED respect unimodality?
%% Construct GED covariance matrices with zeroed crossmodal coefficients

orig_covS = GED.covS;
zeroed_covS = GED.covS;
zeroed_covS(57:end, 1:56) = 0;
zeroed_covS(1:56, 57:end) = 0;

orig_covR = GED.covR;
zeroed_covR = GED.covR;
zeroed_covR(57:end, 1:56) = 0;
zeroed_covR(1:56, 57:end) = 0;

% Note: covR was already regularized, no need to re-do

% Plot as sanity check
% figure; imagesc(zeroed_covS);
% figure; imagesc(zeroed_covR);

%% Perform GED using zeroed covariance matrix

[z_evecs, z_evals] = eig(zeroed_covS, zeroed_covR);
[z_evals, sidx] = sort(diag(z_evals), 'descend');
z_evecs = z_evecs(:, sidx);

% Norm eigenvectors to unit length
for v = 1:size(z_evecs,2)
    z_evecs(:,v) = z_evecs(:,v)/norm(z_evecs(:,v));
end

%% Compute modality dominance values for z_evecs

% Compute modality dominance value for GED eigenvectors (with
% repeating eigenvalues removed)
EEG.nbchan = 56;
for i = 1:size(z_evecs,2)
    e = rms(z_evecs(1:EEG.nbchan,i));
    m = rms(z_evecs(EEG.nbchan+1:end,i));
    md_uni(i) = ((e-m)/(e+m))^2; % this is the magnitude of the dominance ( = unimodality)
end

disp(md_uni); % expectation: all 1 ( = all perfectly unimodal)