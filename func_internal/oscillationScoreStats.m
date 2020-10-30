function [Oscore_rp,fosc_rp,signrep,trendfit] = oscillationScoreStats(cfg, signal)
%{ 
----- oscillationScoreStats -----

DESCRIPTION:
Generates reference signals, either by shuffling or based on trend-fit of 
the original data,and computes O-scores for these reference signals.

INPUTS:
- cfg: struct of settings
- signal: continuous data trace

OUTPUTS:
- Oscore_rp: Oscores of reference signals
- fosc_rp: frequencies of references signals
- signrep: reference signals
- trendfit: trend that was fitted to original data

CONFIG:
- nrep: number of repetitions (required)
- fs: sampling rate of signal (required)
- flim: frequency range of interest (required)
- fpeak: peak frequency of original signal (required)
- keep_trend: switch that determines whether trend is fitted to original
data, if false, data is shuffled (boolean, default = false)
- trend_dist: type of distribution to fit to original data (required for
keep_trend)
- trend_ddt: time step for generation of trend (default = double of fs)
- trend_alpha: significance level of goodness-of-fit test of trend fit
(default = 0.05)
- warnings: verbose switch ('on' or 'off', default = 'on')
- plot: switch for plotting 1st reference trace (boolean; default = false)


Marije ter Wal - 2020
m.j.terwal@bham.ac.uk

%}

%% check config

assert(cfg.nrep>0, 'Invalid number of repetitions')
assert(cfg.fs>0, 'Invalid sampling rate')
assert(isvector(cfg.flim) && length(cfg.flim) == 2, 'Invalid frequency range of interest')
assert(isfield(cfg,'fpeak'), 'No peak frequency specified')
if isnan(cfg.fpeak)
    warning('Peak frequency is NaN, shuffling skipped')
    Oscore_rp = NaN;
    fosc_rp = NaN;
    signrep = [];
    trendfit = [];
    return
    
else
    assert(cfg.fpeak>=cfg.flim(1) && cfg.fpeak<= cfg.flim(2), 'Invalid peak frequency of interest')
end

if ~isfield(cfg,'keep_trend');          cfg.keep_trend = false; end
if cfg.keep_trend 
    if ~isfield(cfg,'trend_dist');      error('No trend distribution was specified'); end
    if ~isfield(cfg,'trend_ddt');       cfg.trend_ddt = 0.5/cfg.fs; end
    if ~isfield(cfg,'trend_alpha');     cfg.trend_alpha = 0.05; end
end

if ~isfield(cfg,'warnings');            cfg.warnings = 'on'; end
if ~isfield(cfg,'plot');                cfg.plot = false; end

%% check signal
if sum(signal)<3
    if strcmp(cfg.warnings, 'on')
        warning('Oscillation score stats could not be computed')
    end
    Oscore_rp = NaN;
    fosc_rp = NaN;
    signrep = [];
    trendfit = [];
    return
end

%% fit distribution

% truncate trace
tstart = max(1,find(signal,1,'first')-1);
tend = min(length(signal), find(signal,1,'last')+1);
signq = signal(tstart:tend);

% check for memory
if length(signal)>30000
    signrep = [];
else
    signrep = zeros(cfg.nrep,length(signal));
end

if isfield(cfg,'keep_trend') && cfg.keep_trend % shuffle while keeping global structure
    
    Dists = cfg.trend_dist;

    pst = ones(length(Dists),1);
    for d = 1:length(Dists)
        
        % fit distribution
        pd = fitdist(find(signq), Dists{d});        
        
        % test fit
        [~,p,stats] = chi2gof(find(signq),'CDF',pd);
        pst(d) = p;
    end
    
    % quality checks
    if isnan(pst) == length(Dists) % all p-vals are NaN
        if strcmp(cfg.warnings, 'on')
            warning('Cannot fit distribution, shuffling data instead.')
        end
        cfg.keep_trend = false;
        
        trendfit = [];
        trendfit.distribution = 'Shuffle';
        trendfit.pd = [];
        trendfit.pval = NaN;
        trendfit.stats = [];
        trendfit.trace = [];
    else
        
        % test the fits
        ids = find(pst>=cfg.trend_alpha); 
        % the null hypothesis for the GOF test is that the data come from
        % the fitted distribution, so we would prefer this not to be
        % rejected
        if isempty(ids) && strcmp(cfg.warnings, 'on')
            warning('No suitable reference distribution found, using the best option.')
        end
        
        % determine the winning distribution
        if length(Dists)>1
            [~,id] = max(pst); 
        else
            id = 1;
        end
        
        % make a trace of the winning trace
        pd = fitdist(find(signq), Dists{id});
        [~,p,stats] = chi2gof(find(signq),'CDF',pd);
        env =  pdf(pd,1:length(signq));
        
        % store fit info
        trendfit = [];
        trendfit.distribution = Dists{id};
        trendfit.pd = pd;
        trendfit.pval = p;
        trendfit.stats = stats;
        trendfit.trace = zeros(size(signal));
        trendfit.trace(tstart:tend) = env/sum(env);
    end
else
    trendfit = [];
    trendfit.distribution = 'Shuffle';
    trendfit.pd = [];
    trendfit.pval = NaN;
    trendfit.stats = [];
    trendfit.trace = [];
end

%% generate reference data and compute Oscore

Oscore_rp = nan(cfg.nrep,1);
fosc_rp = nan(cfg.nrep,1);
signrep = zeros(cfg.nrep,length(signal));

for rp = 1:cfg.nrep
    
    if strcmp(trendfit.distribution, 'Shuffle')
        % create reference traces by shuffling time stamps 
        signrp = zeros(size(signq));
        for i = find(signq)'
            if isnan(cfg.fpeak)
                wind = (cfg.fs/cfg.flim(1));
            else
                wind = (cfg.fs/cfg.fpeak);
            end
           j = max(1, min(length(signrp),round(rand*wind-wind/2) + i )); 
           signrp(j) = signrp(j)+1;
        end
        
    else
        timeline = (tstart)/cfg.fs:cfg.trend_ddt:(tend)/cfg.fs;
        fact = cfg.fs*cfg.trend_ddt;
        envfull = pdf(pd,1:fact:length(signq));
        
        % generate button press time stamps based on trend
        signrpfull = zeros(length(timeline),1);
        signrpfull(rand(size(envfull)) < sum(signal)*envfull/sum(envfull)) = 1;
        
        % reduce the time step
        cfg_mct = [];
        cfg_mct.sd_smooth = [];
        cfg_mct.dt = 1/cfg.fs;
        cfg_mct.removeVal = 0;
        cfg_mct.warnings = 'off';
        
        [signd, timed] = makeContinuousTrace(cfg_mct, timeline(find(signrpfull)));
        
        signrp = zeros(length(signq),1);
        if length(signd)>length(signrp)
            signrp = signd(length(signd)-length(signrp)+1:end);
        else
            signrp(1:length(signd)) = signd;
        end        
          
    end
    
    signdum = zeros(size(signal));
    signdum(tstart:tend) = signrp;
    if ~isempty(signrep)
        signrep(rp,:) = signdum; % store the reference traces
    end
    
    % compute Oscore for reference trace
    [Oscore_rp(rp), fosc_rp(rp),~,~,~] = oscillationScore(cfg, signdum);
    
     if cfg.plot
            figure('Position',[50,50,600,800]);
            subplot(4,1,1); hold on;
            plot([1:length(signq)]/cfg.fs,signq,'color',[0.1 0.5 0.2])
            if isfield(cfg,'keep_trend') && cfg.keep_trend 
            plot([1:length(signq)]/cfg.fs,50*(sum(signq)*env)/sum(env),'color',[0.1 0.5 0.2], 'linewidth',2)
            end
            title('Original data and fit')
            subplot(4,1,2)
            plot([1:length(signq)]/cfg.fs,signrp, 'color', [0.4,0.4,0.4]);
            title('Shuffled data')
            subplot(4,1,3:4); hold on
            plot([1:length(signq)]/cfg.fs,cumsum(signq), 'color',[0.1 0.5 0.2], 'linewidth',2);
            plot([1:length(signq)]/cfg.fs,cumsum(signrp), 'color',[0.4,0.4,0.4], 'linewidth',2);
            title('Cumulative distributions')
            xlabel('Time (s)')            
            
            cfg.plot = false; %we do not want this plot 500 times!
    end
end

end