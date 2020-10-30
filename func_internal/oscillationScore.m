function [oscore,fosc,flim,flimfft,freqs] = oscillationScore(cfg,signal)

%{
----- oscillationScore -----
computes the correlation score as published in:
Muresan et al., J Neurophysiol, 2008
for binary data (spike trains, reaction times)
and continuous data (LFPs, ...)

INPUTS
1. cfg: see below
2. signal: time series with sampling rate fs

CONFIG
---- cfg must contain
- fs in s^-1
- flim: [fmin, fmax] in Hz

---- cfg can also specify
- quantile of data to use (only for spikes/RTs!)
     quantlim: [quantile_min quantile_max] values between 0 and 1 (default: [])
- boolean to determine whether ACH has to be smoothed (might not be necessary for continuous signals)
     smoothach: true or false (default: true)
- width of fast smoothing kernel 
     smoothwind: s (default: determined based on flim)
- width of fast smoothing kernel 
      peakwind: s (default: determined based on flim)
- angular threshold for peak detection
      thresangle in degrees: double (default: 10)
- minimum width of frequency band: 
      minfreqbandwidth: double (default: no minimum)
- frequency of interest: 
      fpeak: Hz (default: [])
- plot: true or false (default: false)
- warnings: 'on' or 'off' (default: on)

OUTPUTS
- oscore: Oscillation score
- fosc: Frequency corresponding to Oscillation score
- flim: Frequency band used for the analysis
- flimfft: Spectrum within the frequency bounds
- freqs: frequency axis corresponding to flimfft

Marije ter Wal - 2020
mj.terwal@bham.ac.uk
%}

%% check config

assert(cfg.fs>0, 'Invalid smapling rate')
assert(isvector(cfg.flim) && length(cfg.flim) == 2, 'Invalid frequency range of interest')

if ~isfield(cfg,'mincycles'); cfg.mincycles = 3; end
if ~isfield(cfg,'thresangle'); cfg.thresangle = 10; end

if ~isfield(cfg,'smoothach'); cfg.smoothach = true; end
if ~isfield(cfg,'smoothwind'); cfg.smoothwind = []; end
if ~isfield(cfg,'peakwind'); cfg.peakwind = []; end

if ~isfield(cfg,'quantlim'); cfg.quantlim = []; end
if ~isfield(cfg,'minfreqbandwidth'); cfg.minfreqbandwidth = []; end
if ~isfield(cfg,'fpeak'); cfg.fpeak = []; end

if ~isfield(cfg,'plot'); cfg.plot = false; end
if ~isfield(cfg,'warnings'); cfg.warnings = 'on'; end

%% check signal 

% checking signal quality
if sum(signal)<3
    if strcmp(cfg.warnings, 'on')
        warning('Dataset is (nearly) empty: Oscillation score could not be computed\n')
    end
    oscore = NaN;
    fosc = NaN;
    flim = [];
    flimfft = [];
    freqs = [];
    return
end

% exclude 'outliers' in discreet data
if ~isempty(cfg.quantlim)
    q = int32(quantile(find(signal),cfg.quantlim));
    signq = signal(max(1,q(1)-1):min(q(2)+1,length(signal)));
else
    % make sure the data do not contain long empty stretches at beginning and end
    tstart = max(1,find(signal,1, 'first')-1); 
    tend = min(length(signal), find(signal,1, 'last')+1);
    signq = signal(tstart:tend);
end

%% check fmin and fmax

fmin = max(cfg.flim(1),cfg.mincycles*(cfg.fs)/(length(signq)));
fmax = min(cfg.flim(2),sum(signq)/(length(signq)/cfg.fs));
flim = [fmin,fmax];

% find the desired with of the ACH: w
w = 2^(1+floor(max(log2(2*cfg.mincycles*cfg.fs/cfg.flim(1)),log2(cfg.fs/2))));

% find s_fast and s_slow
if isempty(cfg.smoothwind)
    cfg.smoothwind = min(0.002,134/(1.5*cfg.flim(2)) /1000); % s - following Muresan et al., 2008
end
if isempty(cfg.peakwind)
    cfg.peakwind = 2*134/(1.5*cfg.flim(1)) /1000; % s
end

if ~isempty(cfg.minfreqbandwidth) && (fmax-fmin) < cfg.minfreqbandwidth
    if strcmp(cfg.warnings, 'on')
        warning('Data not sufficient to meet the minimal frequency bandwidth.')
    end
    oscore = NaN;
    fosc = NaN;
    flimfft = [];
    freqs = [];
    return
end

%% compute oscillation score

% STEP 1: Auto Correlation
ach =  xcorr(signq);

% STEP 2: smooth AC
if cfg.smoothach
    sdwind = round(cfg.smoothwind*cfg.fs);
    % gaussian
    gt = -4*sdwind:4*sdwind;
    wind = 1/(sdwind*sqrt(2*pi)) * exp(-1*(gt.^2 / (2*sdwind^2)));
    ach_smooth = conv(ach,wind,'same');
else
    ach_smooth = ach;
end

% STEP 3: remove peak
sdwind = round(cfg.peakwind*cfg.fs);
gt = -4*sdwind:4*sdwind;
wind = 1/(sdwind*sqrt(2*pi)) * exp(-1*(gt.^2 / (2*sdwind^2)));
ach_slow = conv(ach,wind,'same');
thres = tan(pi*(cfg.thresangle)/180);
scalfact = (length(ach)-1) / ach_slow(ceil(length(ach)/2)) ;

if find(size(ach_slow)>2) == 1 % column array
    diff_ach = flipud(diff(ach_slow(1:ceil(length(ach)/2))));
elseif find(size(ach_slow)>2) == 2 % row array
    diff_ach = fliplr(diff(ach_slow(1:ceil(length(ach)/2))));
else
    error('Signal has the wrong size, consider squeezing!')
end

if isempty(diff_ach)
    if strcmp(cfg.warnings, 'on')
        warning('Dataset is (nearly) empty: Oscillation score could not be computed.')
    end
    oscore = NaN;
    fosc = NaN;
    flimfft = [];
    freqs = [];
    return
end

% find the edge of the ACH peak
peakstart = find( ...
    scalfact*diff_ach<=thres & ... %
    cat(find(size(diff_ach)>1),diff(diff_ach),0)<0, 1, 'first');
if isempty(peakstart)
    if strcmp(cfg.warnings, 'on')
    	warning('Central peak could not be determined.\n')
    end
    phw = 1;
else
    phw = peakstart;
end

% STEP 4: fft and spectral peak
% take the positive lags up to w
tmp = ach_smooth(ceil(length(ach)/2)+phw:end);
if length(tmp) < w
    ach_nopeak = zeros(w,1);
    ach_nopeak(1:length(tmp)) = tmp;
else
    ach_nopeak = tmp(1:w);
end

% get fft and peak within freq limit
cfg.flim = flim;
[peakfreqx, freqs, flimfft] = spectralPeak(cfg, ach_nopeak);

% use a specific frequency of interest
if ~isempty(cfg.fpeak)
    [~,id] = min(abs(freqs-cfg.fpeak));
    peakfreqx = freqs(id);
end

% get fft of whole freq range
cfg.flim = [0,cfg.fs/2];
[~, freqstot, totfft] = spectralPeak(cfg, ach_nopeak);

% STEP 5: compute the score
if isempty(peakfreqx) || isnan(peakfreqx)
    oscore = NaN;
    fosc = NaN;
    flimfft = [];
    freqs = [];
else
    oscore = flimfft(freqs==peakfreqx) / nanmean(totfft);
    fosc = peakfreqx;
end

%% plot
if cfg.plot
    figure('Position',[50,50,600,800]);
    subplot(5,1,1); hold on
    plot([1:length(signal)]/cfg.fs, signal, 'color',[0.5 0.8 0.4]);
    if exist('q','var')
        plot((double(q(1))-1+[1:length(signq)])/cfg.fs, signq, 'color',[0.1 0.5 0.2]);
        plot([double(q(1)),double(q(1))]/cfg.fs,[0,3],':k')
        plot([double(q(2)),double(q(2))]/cfg.fs,[0,3],':k')
    end
    title('Raw signal');
    xlim([0,length(signal)/cfg.fs]);
    ylim([0,3])
    xlabel('Relative time (s)')
    subplot(5,1,2:3); hold on;
    plot([-1*(length(signq)-1):(length(signq)-1)]/cfg.fs,ach,'color',[0.5 0.8 0.4]);
    achsm_nopeak = ach_slow;
    achsm_nopeak(ceil(length(ach)/2)-1-phw: ceil(length(ach)/2)-1+phw) = NaN;
    plot([-1*(length(signq)-1):(length(signq)-1)]/cfg.fs, achsm_nopeak,'color', [0.1 0.5 0.2], 'linewidth',2)
    ylim([0,15])
    xlim([-1*length(signq),length(signq)]/cfg.fs)
    title('ACH')
    xlabel('Lag (s)');
    subplot(5,1,4:5);
    semilogx(freqstot, totfft, 'color', [27,117,187]/255+0.2, 'linewidth',1)
    hold on
    semilogx([freqstot(1),freqstot(end)],[nanmean(totfft),nanmean(totfft)],...
        ':','color', [27,117,187]/255+0.2, 'linewidth',2)
    semilogx(freqs, flimfft, 'color', [27,117,187]/255, 'linewidth',2)
    id = find(ismember(freqs,fosc));
    semilogx(freqs(id),flimfft(id),'o','color',[27,117,187]/255, ...
        'markerfacecolor',[27,117,187]/255, 'linewidth',2)
    plot([flim(1),flim(1)],[0,0.2],':k')
    plot([flim(2),flim(2)],[0,0.2],':k')
    xlim([0.25,cfg.fs/2]);
    
    xlabel('Frequency (Hz)')
    title('Fourier transform of ACH')
    
end

end