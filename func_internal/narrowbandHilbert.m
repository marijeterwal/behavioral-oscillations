function [dathilb,datfilt] = narrowbandHilbert(cfg, data)
%{ 
----- narrowbandHilbert -----

DESCRIPTION:
Filters data in a narrow frequency band and then computes the Hilbert
transform

INPUTS:
- cfg: struct with settings
- data: time series

OUTPUTS:
- dathilb: filtered and hilbert transformed data
- datfilt: filtered data

CONFIG:
- fs: sample rate of data (required)
- freqlim: frequency bounds to use for filtering ([fmin fmax]) (required)
- tol: tolerance for detection of filter stabily
- filtorder: filter order
- demean: switch for demeaning (true or false)
- plot: switch for plotting original, filtered and transformed data (true or false) 

Marije ter Wal - 2020
m.j.terwal@bham.ac.uk

%}


% check config
if ~isfield(cfg, 'fs') || isempty(cfg.fs); error('cfg.fs required'); end
if ~isfield(cfg, 'freqlim') || isempty(cfg.freqlim); error('cfg.freqlim required'); end
if ~isfield(cfg, 'tol') || isempty(cfg.tol); cfg.tol = 1e2; end
if ~isfield(cfg, 'filtorder') || isempty(cfg.filtorder); cfg.filtorder = 2; end
if ~isfield(cfg, 'demean') || isempty(cfg.demean); cfg.demean = true; end
if ~isfield(cfg, 'plot') || isempty(cfg.plot); cfg.plot = false; end

% demean data
if cfg.demean
    dat = data - mean(data);
else
    dat = data;
end

% make filter
if cfg.freqlim(1)<=0
    [b,a] = butter(cfg.filtorder*2,(2*cfg.freqlim(2))/cfg.fs, 'low');
else
    [b,a] = butter(cfg.filtorder,2*[cfg.freqlim]/cfg.fs, 'bandpass');
end

% plot filter
if cfg.plot
    fvtool(b,a, 'Fs',cfg.fs)
end
%
% apply filter to data
datfilt = filtfilt(b,a,dat);

if sum(abs(datfilt)>cfg.tol)
    dathilb = NaN(size(datfilt));
    warning('Filter is unstable')
    return
end

% run hilbert
dathilb = hilbert(datfilt-nanmean(datfilt));

% plot results
if cfg.plot
    figure; plot(dat)
    hold on;
    plot(datfilt,'r')
    plot(angle(dathilb))
    plot(abs(dathilb));
end

end