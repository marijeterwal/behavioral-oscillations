function [signal, tspan] = makeContinuousTrace(cfg,data)
%{ 
----- makeContinuousTrace -----

DESCRIPTION:
Take the descrete events at the time points in data and converts them to a 
(smoothed) continues time series.

INPUTS:
- cfg: struct with settings (see below)
- data: vectore with time points

OUTPUTS:
- signal: time series
- tspan: time axis corresponding to signal

CONFIG:
- dt: time step for continuous trace (required)
- sd_smooth: standard deviation for gaussian smoohting kernel (optional)
- width_block: width of block smoothing kernal (optional)
- removeVal: values to be removed from the data (optional)
- quantlim: extremes to be removed from the data (fractions, optional)
([qmin, qmax])
- warnings: switch for warnings (true or false)


Marije ter Wal - 2020
m.j.terwal@bham.ac.uk

%}


if ~isfield(cfg,'dt'); error('Time step is required'); end
if ~isfield(cfg,'warnings');    cfg.warnings = 'on'; end
   
if isfield(cfg, 'sd_smooth') % gaussian smoothing
    if isempty(cfg.sd_smooth)
        wind = zeros(3,1);
        wind(2) = 1;
    else
        sdwind = round(cfg.sd_smooth / cfg.dt);
        gt = 1:8*sdwind;
        wind = normpdf(gt,4*sdwind, sdwind);
    end
elseif isfield(cfg, 'width_block')
    wind = ones(round(cfg.width_block/cfg.dt),1);
else
    wind = zeros(3,1);
    wind(2) = 1;
end

tspan = 0:cfg.dt:(max(data)+cfg.dt);

% make time series of responses
dum = zeros(length(tspan),1);
for r = 1:length(data)
    if isnan(data(r)) || (isfield(cfg, 'removeVal') && abs(data(r)-cfg.removeVal)<=cfg.dt)
        continue
    end
        
    id = [];
    id = find(int32(tspan*(1/cfg.dt))==int32(round(data(r)*(1/cfg.dt))),1, 'first');
    if isempty(id)
        fprintf(num2str(r))
    end
    dum(id) = dum(id)+ 1;
end

if isfield(cfg, 'quantlim') && ~isempty(cfg.quantlim)
    q = int32(quantile(find(dum),cfg.quantlim));
    dum = dum(q(1):q(2));
    tspan = tspan(q(1):q(2));
end

if sum(dum) ~= length(data) && strcmp(cfg.warnings,'on')
    warning('Datapoints missing'); 
end
signal = conv(dum,wind,'same');
    

end