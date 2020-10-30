function [freq, fxx, corfft] = spectralPeak(cfg, data)

%{ 
----- spectralPeak -----

DESCRIPTION:
Computes the spectrum using the Fast Fourier Transform and finds its peak

INPUTS:
- cfg: struct with settings
- data: time series

OUTPUTS:
- freq: frequency belonging to the spectral peak
- fxx: frequency axis
- corfft: spectrum matching fxx

CONFIG:
- fs: sampling rate (required)
- flim: frequency bounds within which the peak should be found ([fmin,
fmax])
- fcor: switch to use log fit to remove 1/f trend (boolean)
- taper: taper that is applied to the data before fft.

Marije ter Wal - 2020
m.j.terwal@bham.ac.uk

%}

%% check config

if ~isfield(cfg, 'fs') || isempty(cfg.fs)
    error('Sampling rate not specified')
else
    fs = cfg.fs;
end

if ~isfield(cfg,'flim'); cfg.flim = []; end
if ~isfield(cfg,'fcor'); cfg.fcor = false; end
if ~isfield(cfg,'taper'); cfg.taper = 'none'; end

%% compute the fft on data

L = length(data);
if strcmp(cfg.taper, 'none')
    tmp = fft(data);
else
    eval(['tmp = fft(',cfg.taper,'(L).*data);'])
end

ps = abs(tmp(1:round(L/2))/L); %%
ps(2:end-1) = 2*ps(2:end-1);

fx = linspace(0,fs/2,round(L/2));
if isempty(cfg.flim)
    fxx = fx;
    datfft = ps;
else
    f1 = find(fx>=cfg.flim(1),1,'first');
    if f1<1
        f1 = 1;
    end
    f2 = find(fx>=cfg.flim(2), 1, 'first')-1;
    if isempty(f2)
        f2 = length(fx);
    end
    fxx = fx(f1:f2);
    datfft = ps(f1:f2);
end

%% 1/f correction

if cfg.fcor
    % fit 1/f
    warning('This 1/f fitting is very basic, proceed with caution...')
    p = polyfit(log10(fxx), reshape(log10(datfft),size(fxx)), 1);
    
    % correct for 1/f
    corfft = datfft - reshape(10.^(log10(fxx)*p(1)+p(2)),size(datfft));
else
    corfft = datfft;
end

%% find peak

if length(corfft)<3 % can't determine a peak if data is too short
    freq = NaN;
else
    % find the highest peak
    [~,x] = findpeaks(double(corfft),'npeaks',1,'sortstr','descend');
    freq = fxx(x);
end
end