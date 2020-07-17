%% Estimate optimal Granger causality order to determine causality between component time series
% Requires the Multivariate Granger Causality toolbox, available at
% http://www.sussex.ac.uk/sackler/mvgc/
% Toolbox version used: mvgc_v1.0, last updated on March 27 2014
% 
% Author: Marrit Zuure, January 2019
% Drawing from toolbox demo script by Anil Seth

close all; clear;

%% Set paths
dirs = setpaths();

%% Set data import preliminaries
[sublist, ~, ~] = getICs2remove();

%% Set order estimation parameters
ds_factor = 4;    % factor by which to downsample the input data. Recommended target srate: 100 to 500 Hz.
momax = 15;      % max model order (in downsampled time steps). Higher takes longer to fit.

icregmode = 'LWR';   % information criteria regression mode ('OLS', 'LWR' or empty for default)
crit2use = 'BIC';    % using Akaike's (AIC) or Bayes' (BIC) information criterion for order selection

order_est = 'fulltrial';   % Whether to estimate the order WINDOWED or over FULLTRIAL. FULLTRIAL tends to produce larger orders.

% Window settings for windowed order estimation + for GC calculation
wlength_ms = 200;       % in ms
wstep_ms = 50;          % step between windows in ms
wlim_ms = [0 800];      % in ms. Where to start and end windowing (for order estimation only; VAR estimation uses full trial).

% How to choose the optimal order
method = 'max';     % Options: 'mean', 'median', 'max'.
% Recommended: use max order, if not too high to fit the VAR model in a reasonable amount of time.

%% Set VAR model estimation parameters
regmode = 'LWR';    % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

% Significance testing
alpha = 0.05;
mhtc = 'FDR'; % multiple hypotheses correction

%% Set analysis parameters
baseline = [-500 -100]; % in ms
window = [0 800]; % in ms

%% Initialize variable to collect order estimations
num_wins = round(diff(wlim_ms)/wstep_ms);
orders = nan(num_wins, length(sublist));

%% Loop over subjects, estimating model order for each subject (on all trials)
for subno = 1:length(sublist)
    disp(['Processing subject ' num2str(subno) ' of ' num2str(length(sublist)) '...']);
    
    %% Load GED data and previous analysis results   
    GED_filename = [dirs.results sublist{subno} '_GED.mat'];
    load(GED_filename);
    
    ana_filename = [dirs.results sublist{subno} '_ana.mat'];
    load(ana_filename);
    
    %% Per component, remove all-trial ERP from each individual trial,
    % to minimize influence of nonstationarities
    Y = midf.compts;
    for c = 1:midf.num_comps
        ERP = mean(squeeze(midf.compts(c,:,:)),2);
        Y(c,:,:) = bsxfun(@minus, squeeze(Y(c,:,:)), ERP);
    end
    
    %% Downsample data
    % Downsample data by ds_factor to improve order estimation and VAR fitting
    % Reshape component time series to components x time x trials, as toolbox dictates
    X_times = resample(MEEG.times, 1, ds_factor);
    X = nan(midf.num_comps, length(X_times), MEEG.trials);
    for c = 1:midf.num_comps
        X(c,:,:) = resample(squeeze(Y(c,:,:)), 1, ds_factor);
    end
    
    % Determine new sampling rate
    fs = MEEG.srate / ds_factor;
    
    if strcmpi(order_est, 'windowed')
        %% Loop over windows (within window limits)
        % Calculate window indices
        wlim = dsearchn(X_times(:), wlim_ms');
        wlength = round((wlength_ms / 1000) * fs);
        wstep = round((wstep_ms / 1000) * fs);

        for w = 1:num_wins
            wstart = wlim(1) + (w-1)*wstep;
            wend = wlim(1) + (w-1)*wstep + wlength;
            X_w = X(:,wstart:wend,:);
            
            % Convert window data to z-scores
            X_w = zscore(X_w,[],2);

            % Calculate information criteria up to specified maximum model order
            [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X_w,momax,icregmode,0);

            if strcmpi(crit2use, 'AIC')
                orders(w, subno) = moAIC;
            elseif strcmpi(crit2use, 'BIC')
                orders(w, subno) = moBIC;
            else
                error(['''' crit2use ''' not a valid parameter for crit2use (try AIC or BIC).'])
            end
        end
    elseif strcmpi(order_est, 'fulltrial')
        % Calculate information criteria up to specified maximum model order
        [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode,0);
        % Note: Z-scoring X here breaks initialization. However, the
        % component time series in X are based on Z-scored trials, and
        % tsdata_to_infocrit internally demeans X, so this likely does not
        % matter for nonstationarities.
        
        if strcmpi(crit2use, 'AIC')
            orders(1, subno) = moAIC;
        elseif strcmpi(crit2use, 'BIC')
            orders(1, subno) = moBIC;
        else
            error(['''' crit2use ''' not a valid parameter for crit2use (try AIC or BIC).'])
        end
    end
    
end

%% Outside subjects loop: reconcile all model orders and pick one to use

mean_order = nanmean(nanmean(orders));
median_order = nanmedian(reshape(orders,1,[]));
max_order = max(max(orders));
order_sd = nanstd(reshape(orders,1,[]));

disp('Order estimation statistics:');
disp(['Mean: ' num2str(mean_order) ' (' num2str(mean_order * (1000 / fs)) ' ms)']);
disp(['Median: ' num2str(median_order) ' (' num2str(median_order * (1000 / fs)) ' ms)']);
disp(['Max: ' num2str(max_order) ' (' num2str(max_order * (1000 / fs)) ' ms)']);
disp(['SD: ' num2str(order_sd) ' (' num2str(order_sd * (1000 / fs)) ' ms)']);

if strcmpi(method, 'mean')
    order = round(mean_order);
elseif strcmpi(method, 'median')
    order = median_order;
elseif strcmpi(method, 'max')
    order = max_order;
else
    error(['''' method ''' not a valid parameter for method (try MEAN or MEDIAN or MAX).'])
end

disp(['Choosing ' method ' order of ' num2str(order) ' (' num2str(order * (1000 / fs)) ' ms)']);

%% Second loop over subjects, using model order to estimate VAR model and calculate GC
for subno = 1:length(sublist)
    disp(['Processing subject ' num2str(subno) ' of ' num2str(length(sublist)) '...']);
    
    %% Load GED data and previous analysis results
    GED_filename = [dirs.results sublist{subno} '_GED.mat'];
    load(GED_filename);
    
    ana_filename = [dirs.results sublist{subno} '_ana.mat'];
    load(ana_filename);
    
    clear gc;
    
    %% Per component, remove all-trial ERP from each individual trial, 
    % to minimize influence of nonstationarities
    Y = midf.compts;
    for c = 1:midf.num_comps
        ERP = mean(squeeze(midf.compts(c,:,:)),2);
        Y(c,:,:) = bsxfun(@minus, squeeze(Y(c,:,:)), ERP);
    end
    
    %% Downsample data
    % Downsample data by ds_factor to improve order estimation and VAR fitting
    % Reshape component time series to components x time x trials, as toolbox dictates
    X_times = resample(MEEG.times, 1, ds_factor);
    X = nan(midf.num_comps, length(X_times), MEEG.trials);
    for c = 1:midf.num_comps
        X(c,:,:) = resample(squeeze(Y(c,:,:)), 1, ds_factor);
    end
    
    % Determine new sampling rate
    fs = MEEG.srate / ds_factor;
    
    %% Loop over windows (across full trial)
    % Calculate window indices
    num_wins = round((X_times(end)-X_times(1))/wstep_ms);
    wlength = round((wlength_ms / 1000) * fs);
    wstep = round((wstep_ms / 1000) * fs);
    
    % Initialize matrix to collect pairwise GC data
    gc.F_pw_cC = nan(midf.num_comps, midf.num_comps, num_wins);
    gc.F_driver_cC = nan(midf.num_comps, num_wins);
    gc.F_receiver_cC = nan(midf.num_comps, num_wins);
    gc.F_pw_cI = nan(midf.num_comps, midf.num_comps, num_wins);
    gc.F_driver_cI = nan(midf.num_comps, num_wins);
    gc.F_receiver_cI = nan(midf.num_comps, num_wins);
    
    for w = 1:num_wins
        wstart = 1 + (w-1)*wstep;
        wend = 1 + (w-1)*wstep + wlength;
        X_w = X(:,wstart:wend,:);
        
        % Convert window data to z-scores
        X_w = zscore(X_w,[],2);
        
        X_cC = X_w(:,:,trialtype==1);
        X_cI = X_w(:,:,trialtype==2);
        
        %% Estimate VAR model of selected order from data.
        ptic('\n*** tsdata_to_var... ');
        [A_cC,SIG_cC] = tsdata_to_var(X_cC,order,regmode);
        ptoc;
        
        ptic('\n*** tsdata_to_var... ');
        [A_cI,SIG_cI] = tsdata_to_var(X_cI,order,regmode);
        ptoc;
        
        % Check for failed regression
        assert(~isbad(A_cC),'VAR estimation failed');
        assert(~isbad(A_cI),'VAR estimation failed');
        
        %% Autocovariance calculation
        ptic('*** var_to_autocov... ');
        [G_cC,info_cC] = var_to_autocov(A_cC,SIG_cC,acmaxlags);
        ptoc;
        
        ptic('*** var_to_autocov... ');
        [G_cI,info_cI] = var_to_autocov(A_cI,SIG_cI,acmaxlags);
        ptoc;
        
        var_info(info_cC,true); % report results (and bail out on error)
        var_info(info_cI,true); % report results (and bail out on error)
        
        %% Granger causality calculation: time domain
        
        % Calculate pairwise conditional causalities for this window
        ptic('*** autocov_to_pwgc... ');
        F_pw_cC = autocov_to_pwcgc(G_cC);
        ptoc;
        
        ptic('*** autocov_to_pwgc... ');
        F_pw_cI = autocov_to_pwcgc(G_cI);
        ptoc;
        
        % Significance tests - pairwise
        pval_pw_cC = mvgc_pval(F_pw_cC,order,size(X_w,2),size(X_w,3),1,1,midf.num_comps-2,'F'); % take careful note of arguments!
        sig_pw_cC  = significance(pval_pw_cC,alpha,mhtc);
        pval_pw_cI = mvgc_pval(F_pw_cC,order,size(X_w,2),size(X_w,3),1,1,midf.num_comps-2,'F'); % take careful note of arguments!
        sig_pw_cI  = significance(pval_pw_cI,alpha,mhtc);
        
        F_driver_cC = nan(1,midf.num_comps);
        F_receiver_cC = nan(1,midf.num_comps);
        F_driver_cI = nan(1,midf.num_comps);
        F_receiver_cI = nan(1,midf.num_comps);
        
        % Calculate total driving and receiving strength per component for this window:
        for c = 1:midf.num_comps
            F_driver_cC(c) = autocov_to_mvgc(G_cC,c,[1:c-1 c+1:midf.num_comps]); % total driving strength for this component (same as sum pw)
            F_receiver_cC(c) = autocov_to_mvgc(G_cC,[1:c-1 c+1:midf.num_comps],c); % total receiving strength for this component (same as sum pw)
            
            F_driver_cI(c) = autocov_to_mvgc(G_cI,c,[1:c-1 c+1:midf.num_comps]); % total driving strength for this component (same as sum pw)
            F_receiver_cI(c) = autocov_to_mvgc(G_cI,[1:c-1 c+1:midf.num_comps],c); % total receiving strength for this component (same as sum pw)
            
            % Significance tests
            pval_driver_cC(c) = mvgc_pval(F_driver_cC(c),order,size(X_w,2),size(X_w,3),midf.num_comps-1,1,0,'chi2'); % take careful note of arguments!
            sig_driver_cC(c)  = significance(pval_driver_cC(c),alpha,mhtc);
            pval_receiver_cC(c) = mvgc_pval(F_receiver_cC(c),order,size(X_w,2),size(X_w,3),1,midf.num_comps-1,0,'chi2'); % take careful note of arguments!
            sig_receiver_cC(c)  = significance(pval_receiver_cC(c),alpha,mhtc);
            
            pval_driver_cI(c) = mvgc_pval(F_driver_cI(c),order,size(X_w,2),size(X_w,3),midf.num_comps-1,1,0,'chi2'); % take careful note of arguments!
            sig_driver_cI(c)  = significance(pval_driver_cI(c),alpha,mhtc);
            pval_receiver_cI(c) = mvgc_pval(F_receiver_cI(c),order,size(X_w,2),size(X_w,3),1,midf.num_comps-1,0,'chi2'); % take careful note of arguments!
            sig_receiver_cI(c)  = significance(pval_receiver_cI(c),alpha,mhtc);
        end
         
        % Normalize by num_comps - 1 (number of incoming/outgoing
        % connections per node)
        F_driver_cC = F_driver_cC ./ (midf.num_comps-1);
        F_driver_cI = F_driver_cI ./ (midf.num_comps-1);
        F_receiver_cC = F_receiver_cC ./ (midf.num_comps-1);
        F_receiver_cI = F_receiver_cI ./ (midf.num_comps-1);
        
        % Check for failed GC calculation
        assert(~isbad(F_pw_cC,false),'Pairwise GC calculation failed');
        assert(~isbad(F_pw_cI,false),'Pairwise GC calculation failed');
        
        %% Collect data into variable to save
        gc.F_pw_cC(:,:,w) = F_pw_cC;
        gc.F_driver_cC(:,w) = F_driver_cC;
        gc.F_receiver_cC(:,w) = F_receiver_cC;
        gc.F_pw_cI(:,:,w) = F_pw_cI;
        gc.F_driver_cI(:,w) = F_driver_cI;
        gc.F_receiver_cI(:,w) = F_receiver_cI;
        gc.F_pw(:,:,w) = mean(cat(4, F_pw_cC, F_pw_cI),4); % mean of cC and cI trial models (instead of fitting an independent one)
        gc.F_driver(:,w) = mean(cat(3, F_driver_cC, F_driver_cI),3);
        gc.F_receiver(:,w) = mean(cat(3, F_receiver_cC, F_receiver_cI),3);
        
        gc.sig_pw_cC(:,:,w) = sig_pw_cC;
        gc.pval_pw_cC(:,:,w) = pval_pw_cC;
        gc.sig_pw_cI(:,:,w) = sig_pw_cI;
        gc.pval_pw_cI(:,:,w) = pval_pw_cI;
        gc.sig_driver_cC(:,w) = sig_driver_cC;
        gc.pval_driver_cC(:,w) = pval_driver_cC;
        gc.sig_receiver_cC(:,w) = sig_receiver_cC;
        gc.pval_receiver_cC(:,w) = pval_receiver_cC;
        gc.sig_driver_cI(:,w) = sig_driver_cI;
        gc.pval_driver_cI(:,w) = pval_driver_cI;
        gc.sig_receiver_cI(:,w) = sig_receiver_cI;
        gc.pval_receiver_cI(:,w) = pval_receiver_cI;
        
        % Clear variables
        clear pval_driver_cC pval_driver_cI pval_receiver_cC pval_receiver_cI sig_driver_cC sig_driver_cI sig_receiver_cC sig_receiver_cI;
    end
    
    %% Compute GC baseline, task modulation, and conflict modulation
    
    winmid = [];
    for i = 1:num_wins
        winmid = [winmid; MEEG.times(1) + wlength_ms/2 + (i-1)*wstep_ms];
    end
    
    widx = dsearchn(winmid(:), window');
    bidx = dsearchn(winmid(:), baseline');
    
    gc.driver_baseline = mean(gc.F_driver(:, bidx(1):bidx(2), :),2);
    gc.receiver_baseline = mean(gc.F_receiver(:, bidx(1):bidx(2), :),2);
    
    gc.driver_task = mean(gc.F_driver(:, widx(1):widx(2),:),2) ./ gc.driver_baseline;
    gc.receiver_task = mean(gc.F_receiver(:, widx(1):widx(2),:),2) ./ gc.receiver_baseline;
    % GC task modulation is normalized by an all-trial -500:-100 ms baseline per component, 
    % conceptually similar to how theta power task modulation is computed.
    
    gc.driver_conflict = mean(gc.F_driver_cI(:, widx(1):widx(2),:),2) ./ gc.driver_baseline ...
        - mean(gc.F_driver_cC(:,widx(1):widx(2),:),2) ./ gc.driver_baseline;
    gc.receiver_conflict = mean(gc.F_receiver_cI(:, widx(1):widx(2),:),2) ./ gc.receiver_baseline ...
        - mean(gc.F_receiver_cC(:,widx(1):widx(2),:),2) ./ gc.receiver_baseline;
    
    %% Collect extra variables into struct to save
    % Didn't specify these as gc fields upfront, as gc (if it exists already)
    % can be loaded in the second subjects loop in this file and would have 
    % overwritten the parameters.
    gc.ds_factor = ds_factor;
    gc.momax = momax;
    gc.crit2use = crit2use;
    gc.order_est = order_est;
    gc.wlength = wlength_ms;
    gc.wstep = wstep_ms;
    gc.wlim = wlim_ms;
    gc.method = method;
    gc.baseline = baseline;
    gc.window = window;
    gc.num_wins = num_wins;
    gc.orders = orders;
    gc.mean_order = mean_order;
    gc.median_order = median_order;
    gc.max_order = max_order;
    gc.order_sd = order_sd;
    gc.order = order;
    gc.winmid = winmid;
    
    %% Append data to analysis results (for later plotting)
    
    disp('Saving results to file...');
    ana_filename = [dirs.results sublist{subno} '_ana.mat'];
    save(ana_filename, 'gc', '-append');
end

%%
disp('Run completed successfully.');