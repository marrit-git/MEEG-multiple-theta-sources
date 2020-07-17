function [ tf ] = tfdecomp( data, MEEG, min_freq, max_freq, num_frex, method, varargin )
%TFDECOMP decomposes a matrix of time series data into time-frequency power.
%   TFDECOMP(data, MEEG, min_freq, max_freq, num_frex, method, times2save,
%   baseline, trialtype) uses sliding wavelet convolution to extract 
%   power over time at each of the given frequencies. Pass method 'means' 
%   to return TF decomposition with mean over trials per condition; pass 
%   method 'trials' to get the TF decomposition of individual trials 
%   (without dB scaling, as logarithmic scales are sensitive to 
%   fluctuations); pass method 'phase' to get just the phase angle time series.
%
%   FUNCTION INPUTS
%       data            The data to decompose; sensors x time x trials
%       MEEG            Struct containing at least fields pnts, times,
%                          trials, srate
%       min_freq        Lowest frequency (in Hz) at which to decompose
%       max_freq        Highest frequency (in Hz) at which to decompose
%       num_frex        Number of frequencies at which to decompose
%       method          'means', 'trials', or 'phase'
%
% If METHOD is 'means' or 'trials', also pass the following:
%       times2save      Array of times (in ms) at which to decompose
%       baseline        2x1 vector of baseline start and end time (in ms)
%
% If METHOD is 'means', also pass the following:
%       trialtype       Integer array listing task condition per trial. If
%                       not passed, only mean over all trials is returned.
%
%   FUNCTION OUTPUTS
%       tf (method 'means')  Power time series; sensors x frequencies x 
%                               time x conditions+1. First condition is 
%                               mean over all trials; other conditions are 
%                               means over each condition, sorted by the 
%                               integer value in trialtype.
%       tf (method 'trials') Amplitude time series for individual trials;
%                               sensors x frequencies x time x trials
%       tf (method 'phase')  Phase time series per frequency per trial; 
%                               sensors x frequencies x timepoints x trials
%
%% Pre-check input validity
if ~strcmpi(method, 'means') && ~strcmpi(method, 'trials') && ~strcmpi(method, 'phase')
    error(['Call to tfdecomp failed: ''' method ''' invalid value for input parameter ''method'' (try MEANS or TRIALS or PHASE)']);
end

%% Handle variable number of arguments
try
    times2save = varargin{1};
    baseline = varargin{2};
    trialtype = varargin{3};
end

%% Extract extra variables from inputs
frex = logspace(log10(min_freq),log10(max_freq),num_frex);

% Convert time points from ms to index
if ~strcmpi(method, 'phase')
    tidx = dsearchn(MEEG.times(:),times2save');
    bidx = dsearchn(MEEG.times(:),baseline');
end

% Extract the number of conditions
if strcmpi(method, 'means')
    if ~exist('trialtype', 'var') % if trialtype was not passed
        conditions = [];
    else
        conditions = unique(trialtype, 'sorted');
    end
end

%% Create Morlet wavelets: multiply Gaussians by sine waves

% Gaussian width
if num_frex == 1 && min_freq > 5 && min_freq < 7 % When only decomposing at theta, use a 
    % slightly narrower wavelength width than the one returned by the
    % equation below. Chosen to correspond to theta in logspace(log10(4),
    % log10(10), num_frex) where min_freq is 2, max_freq is 30, num_frex is 40.
    s = 2*( 6.25 ./(2*pi*frex) ).^2;
else
    s = 2*( (logspace(log10(4),log10(10),num_frex))./(2*pi*frex) ).^2;
    % Note that the wavelet width can be different for a given frequency for
    % different lengths of frex, minimum frex, or maximum frex.
end

% Create wavelets
wt = -2:1/MEEG.srate:2;
nWave = length(wt);
halfw = floor(nWave/2)+1;
nConv = MEEG.pnts*MEEG.trials + nWave - 1;

waveX = zeros(num_frex,nConv);

for fi=1:num_frex
    tmp = exp(1i*2*pi*frex(fi)*wt) .* exp( -wt.^2/s(fi) );
    waveX(fi,:) = fft( tmp./max(tmp) ,nConv);
end

%% Get power spectra by wavelet convolution
    
% Initialize time-frequency matrix
if strcmpi(method, 'means')
    tf = zeros(size(data,1),num_frex,length(times2save),length(conditions)+1);
elseif strcmpi(method, 'trials')
    tf = zeros(size(data,1),num_frex,length(times2save),MEEG.trials);
elseif strcmpi(method, 'phase')
    tf = zeros(size(data,1),num_frex,MEEG.pnts,MEEG.trials);
end

% Loop over time series (may be sensors or components, depending on input data)
for tsi = 1:size(data,1)
    disp(['TF decomposing time series ' num2str(tsi) ' of ' num2str(size(data,1)) ' (method ''' method ''')...'])
    
    % Data spectrum from this time series
    dataX = fft( reshape(data(tsi,:,:),1,[]) ,nConv);
    
    % Loop over frequencies
    for fi=1:num_frex
        % Get power time series
        as = ifft( dataX.*waveX(fi,:) );
        as = reshape( as(halfw:end-halfw+1) , [MEEG.pnts MEEG.trials] ).^2;
        
        if strcmpi(method, 'means') % return TF decomposition mean per condition
            as  = abs(as); % extract magnitude from complex analytic signal
            base = mean(mean(as(bidx(1):bidx(2),:),2),1); % baseline
            
            tf(tsi, fi, :, 1) = 10*log10( mean(as(tidx,:),2) / base ); % all trials
            if length(conditions) > 1
                for tti = 1:length(conditions) % trials belonging to a given condition
                    cond = conditions(tti);
                    tf(tsi, fi, :, tti+1) = 10*log10( mean(as(tidx,trialtype==cond),2) / base );
                end
            end
        elseif strcmpi(method, 'trials') % return TF decomposition per trial
            as  = abs(as); % extract magnitude from complex analytic signal
            base = mean(mean(as(bidx(1):bidx(2),:),2),1); % baseline
            
            tf(tsi,fi,:,:) = as(tidx,:) / base;
        elseif strcmpi(method, 'phase') % return phase angle of the complex analytic signal
            tf(tsi,fi,:,:) = angle(as);
        end
    end
end
end



