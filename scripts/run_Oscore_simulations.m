
%{
This script produces simulated behavioral data and assesses the performance 
of the O-score and Z-scoring procedures. 
It produces plots for Figure S5 from:

Ter Wal, M., Linde Domingo, J., Lifanov, J., Roux, F., Kolibius, L., 
Gollwitzer, S., Lang, J., Hamer, H., Rollings, D., Sawlani, V., 
Chelvarajah, R., Staresina, B., Hanslmayr, S., Wimber, M. 
Theta rhythmicity governs the timing of behavioural and hippocampal 
responses in humans specifically during memory- dependent tasks. 
bioRxiv 2020.11.09.374264; DOI: https://doi.org/10.1101/2020.11.09.374264

For details on the O-score method, please refer to:
Muresan, R.C., Jurjut, O.F., Moca, V. V., Singer, W., and Nikolic, D. (2008). 
The Oscillation Score: An Efficient Method for Estimating Oscillation Strength 
in Neuronal Activity. J. Neurophysiol. 99, 1333–1353.


This script depends on the following functions:

Internal - Oscore:
- makeContinuousTrace
- oscillationScore
- oscillationScoreStats
- spectralPeak


Marije ter Wal - 2020
m.j.terwalm@bham.ac.uk

%}


%% settings

% choose the characteristics of the simulated data:
taskphase = 'visual'; % choose from 'visual','encoding' of 'retrieval'
osc_freq = 5; % Hz; frequency of oscillatory modulation of responses
osc_amps = 0:0.1:1; % strength of oscillatory modulation (between 0 and 1)

saveresults = false;
savepath = 'YourPathGoesHere';

% simulation time steps (same as for behavioral data)
ddt = 0.0005; % timestep, seconds
dt = 0.001; % timestep, seconds
seed = 2;

% task phase specific settings - these are based on Table S3
if sum(strcmp(taskphase, {'visual','v'}))
    % for VISUAL
    Ttot = 3;% total distribution time, in seconds
    nresp = 215; % average number of responses per subject
    toffset = 0.4; % minimum reaction time
    sdresp = 54; % SD of number of responses per subject
    nsubj = 95; % number of included subjects
elseif sum(strcmp(taskphase, {'encoding','e'}))
    % for ENCODING
    Ttot = 8;
    nresp = 66;
    toffset = 1;
    sdresp = 34;
    nsubj = 190;
elseif sum(strcmp(taskphase, {'retrieval','r'}))
    % for RETRIEVAL
    Ttot = 8;
    nresp = 151;
    toffset = 0.4;
    sdresp = 54;
    nsubj = 70;
end

% for oscore (same as for behavioral data)
cfgt_oscore = [];
cfgt_oscore.quantlim = [0.05,0.95];
cfgt_oscore.mincycles = 3;
cfgt_oscore.fs = 1/dt;
cfgt_oscore.smoothwind = 0.002;
cfgt_oscore.peakwind = 0.008;
cfgt_oscore.warnings = 'off';
cfgt_oscore.taper = 'hanning';
cfgt_oscore.fcor = false;
cfgt_oscore.flim = [0.5,40];
cfgt_oscore.nrep = 500;
cfgt_oscore.keep_trend = true;
cfgt_oscore.trend_dist = {'Gamma'};
cfgt_oscore.trend_ddt = ddt;
cfgt_oscore.trend_alpha = 0.05; % alpha for goodness of fit
cfgt_oscore.thresangle = 10;
cfgt_oscore.plot = false;

cfg_mct = [];
cfg_mct.sd_smooth = [];
cfg_mct.dt = dt;
cfg_mct.removeVal = 1;
cfg_mct.warnings = 'off';

% for plotting
cl = [];
cl.red = [190 30 45]/255;
cl.blue = [27,117,187]/255;%[56 100 176]/256;
cl.gray = [0.6,0.6,0.6];

%% loop over osc. amplitudes

oscore = nan(length(osc_amps),nsubj);
oscoreZ = nan(length(osc_amps),nsubj);
fosc = nan(length(osc_amps),nsubj);
RTs = cell(length(osc_amps),nsubj);

phases = {};

for a = 1:length(osc_amps)
    
    fprintf('Calculating amp %i of %i:', a,length(osc_amps))
    
    rng(seed)
    osc_amp = osc_amps(a);
    
    for s = 1:nsubj
        
        cfg_oscore = cfgt_oscore;
        
        fprintf(' %i ', s)
        
        % --------- simulate spike trains
        
        % randomize the number of responses
        snresp = max(10,round(randn*sdresp) + nresp);
        % and the length of the trace
        DtotT = (0.5+rand)*Ttot;
        
        timeline = 0:ddt:DtotT;
        
        % create the oscillatory trace
        osctrace = 1+osc_amp*sin(2*pi*osc_freq.*timeline);
        
        % and create the trend trace
        % these are based on Table S3
        if sum(strcmp(taskphase, {'visual','v'}))
            % visual
            trendtrace = gampdf(timeline,1+rand,(1+rand)/4); % gamma
        elseif sum(strcmp(taskphase, {'encoding','e'}))
            % encoding
            trendtrace = normpdf(timeline,1.5+rand,2+rand); % normal
        elseif sum(strcmp(taskphase, {'retrieval','r'}))
            % retrieval
            trendtrace = lognpdf(timeline,rand,1+0.5*rand); % lognormal         
        end
            
        % combine trend and oscillation
        ratetrace = trendtrace .* osctrace * snresp;
        
        % generate the button presses
        ddata = zeros(size(timeline));
        ddata(rand(size(timeline)) < ratetrace*ddt) = 1;
        RTs{a,s} = toffset+timeline(find(ddata));
        
        [data,timelinet] = makeContinuousTrace(cfg_mct,RTs{a,s});

        % --------- calculate o-score
        cfg_oscore.fpeak = [];
        [oscore(a,s),fosc(a,s),~,~,~] = oscillationScore(cfg_oscore,data);
        
        % o-score stats
        cfg_oscore.fpeak = fosc(a,s);
        [oscore_rep, fosc_rep,~,~] = oscillationScoreStats(cfg_oscore, data);
        oscoreZ(a,s) = (log(oscore(a,s)) - nanmean(log(oscore_rep))) ./ nanstd(log(oscore_rep));
        
    end
    fprintf(' Done!\n')
end

%% Figure S5: plot o-scores and frequencies

alpha = 0.05;
thres = norminv(1-alpha,0,1);

% stats
dum = oscoreZ - thres;
vargroup = nansum((dum - repmat(nanmean(dum,2),[1,nsubj])).^2,2)/(nsubj-1);
Tpop_Os = nanmean(dum,2) ./ sqrt(vargroup/(nsubj-1));
pval = tcdf(Tpop_Os,nsubj-1,'upper');
sign = nan(length(osc_amps),1);
sign(pval<alpha/5) = 1; % use same statistical threshold as for the data

figure;
% o-scores
subplot(311); hold on;
plot([-0.05,1.05], [0,0],'k')
plot([-0.05,1.05], [thres,thres],':k')
errorbar(osc_amps, nanmean(oscoreZ,2), nanstd(oscoreZ,0,2),'color',cl.blue)
plot(osc_amps, nanmean(oscoreZ,2),'color',cl.blue,'linewidth',2)
plot(osc_amps, 0.5*sign, '*','color',cl.gray)
ylim([0,4])
xlim([-0.05,1.05])
ylabel('O-score (Z-scored)')

% frequencies
subplot(312); hold on;
dum = round(fosc);
dum(oscoreZ<thres) = NaN;
cmap = brewermap(100,'Blues');
hb = nan(length(osc_amps),21);
for o = 1:length(osc_amps)
   hb(o,:) = hist(dum(o,:),0:20);
end
imagesc(osc_amps,0:20,hb'); colormap(cmap); axis xy;
caxis([0,size(dum,2)/4]); % color bar between 0 and 25% of subjects
plot([-0.05,1.05],[osc_freq,osc_freq]-1,':','color',cl.red,'linewidth',1)
plot([-0.05,1.05],[osc_freq,osc_freq]+1,':','color',cl.red,'linewidth',1)
ylabel('Peak frequency (Hz)')
xlim([-0.05,1.05])
ylim([-0.5,20.5])

% fraction significant
subplot(313); hold on;
bar(osc_amps, sum(oscoreZ>=thres,2)/nsubj, 0.2, 'facecolor',cl.blue, 'edgecolor','none')%,'linewidth',2, 
ylabel('Fraction significant')
xlim([-0.05,1.05])
ylim([0,1])

xlabel('Amplitude of oscillation (%)')

if saveresults
    % save figures
    saveas(gcf, [savepath, 'sim_',taskphase, '_',num2str(osc_freq),'Hz',...
        '_mincyc',num2str(cfgt_oscore.mincycles),...
        '_smoothwind', num2str(cfgt_oscore.smoothwind*1000),...
        '_nrep',num2str(cfgt_oscore.nrep),'_',cfg_oscore.taper,'.fig'])
end

%% save data

if saveresults
    % save data
    save([savepath, 'sim_',taskphase, '_',num2str(osc_freq),'Hz',...
        '_mincyc',num2str(cfgt_oscore.mincycles),...
        '_smoothwind', num2str(cfgt_oscore.smoothwind*1000),...
        '_nrep',num2str(cfgt_oscore.nrep),'_',cfg_oscore.taper, '.mat'],...
        'osc_amps','osc_freq','nsubj',...
        'oscoreZ','fosc')
end
