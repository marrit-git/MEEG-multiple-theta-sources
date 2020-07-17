function [ perm_evals, perm_settings, perm_evecs, perm_covS ] = permuteGED( MEEG, covS_window, covR_window, kernel, num_perm, method)
%PERMUTEGED Generates pool of permutation-driven eigenvalues for GED
%significance thresholding.
%   PERMUTEGED randomly shuffles each (filtered, unfiltered) trial data
%   into (covS, covR) and applies a generalized eigendecomposition to the 
%   resulting covariance matrices. This is repeated NUMPERM times. The
%   resulting eigenvalues can be used for Z-scoring and significance
%   thresholding. 
%   Function can return either the complete pool of eigenvalues, or the 
%   pool containing only the highest eigenvalue on each iteration. The
%   latter can be used to control multiple comparisons correction
%   using the MinP/MaxT method.
%
%   FUNCTION INPUTS
%       MEEG            The data matrix; sensors x timepoints x trials
%       covS_window     Time window (1x2) to use for constructing S, in ms
%       covR_window     Time window (1x2) to use for constructing R, in ms
%       filter_lo       Filter lower bound, in Hz
%       filter_hi       Filter upper bound, in Hz
%       filter_order    Filter order
%       num_perm         The number of permutations to complete
%       method          'ALL', 'MAX', or '95', to return complete pool of
%                           eigenvalues, only the highest from each
%                           iteration, or the 95th percentile for that
%                           specific eigenvalue.
%
%   FUNCTION OUTPUTS
%       perm_evals      Array of eigenvalues; either all (num_perm*sensors x
%                           1) or max (num_perm x 1)
%       perm_settings   List containing filter settings and window, for 
%                       later comparison to GED settings
%       perm_evecs      num_perm x sensors x sensors matrix of eigenvectors
%       perm_covS       num_perm x sensors x sensors covariance matrix

%% Pre-check input validity
if ~strcmpi(method, 'all') && ~strcmpi(method, 'max') && ~strcmpi(method, '95')
    error(['Call to permuteGED failed: ''' method ''' invalid value for input parameter ''method'' (try ''ALL'' or ''MAX'' or ''95'')']);
end

%% Save some input settings, for future comparison to GED settings
perm_settings = [covS_window, covR_window, kernel];

%% Filter data
disp('Filtering data...');
fft_len = size(MEEG.data,2)+length(kernel)-1; % number of time points/frequencies for fft to return
trim_len = (length(kernel)-1)/2; % number of time points to trim from start and end of result
fdata = 2*real( ifft( bsxfun(@times,fft(MEEG.data,fft_len,2),fft(kernel, fft_len)) ,[],2) );
fdata = reshape(fdata(:,trim_len+1:end-trim_len,:),size(MEEG.data,1), MEEG.pnts, MEEG.trials);

%% Convert time windows for GEDs from ms to timepoint index
tidx_S = dsearchn(MEEG.times(:),covS_window');
tidx_R = dsearchn(MEEG.times(:),covR_window');

%% Loop over permutations

perm_evals = [];
perm_evecs = nan(num_perm, MEEG.nbchan, MEEG.nbchan);
perm_covS = nan(num_perm, MEEG.nbchan, MEEG.nbchan);

for permi = 1:num_perm
    disp(['Permutation testing GED ' num2str(permi) ' of ' num2str(num_perm) '...']);
    
    covSperm = zeros(MEEG.nbchan);
    covRperm = zeros(MEEG.nbchan);
    
    for triali=1:MEEG.trials
        Sdata = fdata(:,:,triali);
        Sdata = zscore(Sdata(:,tidx_S(1):tidx_S(2)));
        
        Rdata = MEEG.data(:,:,triali);
        Rdata = zscore(Rdata(:,tidx_R(1):tidx_R(2)));
        
        if rand<.5
            covSperm = covSperm + Sdata*Sdata' / size(Sdata,1)';
            covRperm = covRperm + Rdata*Rdata' / size(Rdata,1) ;
        else
            % randomly swap which data go into which covariance matrix
            covSperm = covSperm + Rdata*Rdata' / size(Rdata,1);
            covRperm = covRperm + Sdata*Sdata' / size(Sdata,1);
        end
    end
    
    covSperm = covSperm ./ triali;
    covRperm = covRperm ./ triali;
    
    % Regularization
    g = 0.01;
    covRrperm = (1-g)*covRperm + g*mean(eig(covRperm))*eye(size(MEEG.data,1));
    
    % Do GED on S_perm vs R_perm
    [evecs, evals] = eig(covSperm, covRrperm);
    evals = sort(diag(evals), 'descend');
    if strcmpi(method, 'max')
        perm_evals = [perm_evals; evals(1)]; % add maximum eigenvalue to pool ("max-T method")
    elseif strcmpi(method, 'all')
        perm_evals = [perm_evals; evals]; % add all eigenvalues to pool
    elseif strcmpi(method, '95')
        perm_evals = [perm_evals, evals]; % add all eigenvalues to pool; extract only 95th percentile after running all permutations
    end
    perm_covS(permi,:,:) = covSperm;
    perm_evecs(permi,:,:) = evecs;
end
if strcmpi(method, '95')
    temp_evals = zeros(1,length(covSperm));
    for i = 1:length(covSperm)
        temp_evals(i) = prctile(perm_evals(i,:), 95); % for each eigenvalue, extract 95th percentile from full set of permutations
    end
    perm_evals = temp_evals;
end

end

