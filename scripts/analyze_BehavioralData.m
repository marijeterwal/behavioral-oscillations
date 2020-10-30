%{
This script analyses the behavioral data and produces plots for Figures 2-4 from:

Ter Wal, M., Linde Domingo, J., Lifanov, J., Roux, F., Kolibius, L., Gollwitzer, S., 
Lang, J., Hamer, H., Rollings, D., Sawlani, V., Chelvarajah, R., Staresina, B., 
Hanslmayr, S., Wimber, M., (in preparation), Theta rhythmicity governs the 
timing of behavioural and hippocampal responses in humans specifically during 
memory- dependent tasks.

The data are available here:
Ter Wal, M., Linde Domingo, J., Lifanov, J., Roux, F., Kolibius, L., 
Gollwitzer, S., Lang, J., Hamer, H., Rollings, D., Sawlani, V., Chelvarajah, 
R., Staresina, B., Hanslmayr, S., Wimber, M., (2020), Data for: Theta 
rhythmicity governs the timing of behavioural and hippocampal responses in 
humans specifically during memory- dependent tasks. figshare. Collection. 
DOI: 10.6084/m9.figshare.c.5192567


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
- narrowbandHilbert

Internal - plotting
- simplebox

External:
- CircStat Toolbox 2012a: circ_vtest. 
Please refer to: Berens, P. (2009). CircStat: A MATLAB Toolbox for Circular 
Statistics. J. Stat. Softw. 31.
This function can be found here: https://github.com/circstat/circstat-matlab
- rotateXLabels


Marije ter Wal - 2020
m.j.terwalm@bham.ac.uk

%}

%% settings

% select task phase:
% choose from: 'encoding', 'ret&catch', 'retrieval', 'catch', 'visual'
taskphase = 'retrieval';

datapath = 'Z:\Marije\Data_RT_logs\Logs_ForPublication\';
storepath = 'Z:\Marije\Project_RT\';

storename = '_test';

% set random number generator for reproducibility
rng(1)

% save?
save_figures = false; %true;
save_data = false; %true;

% run phase analyses?
analysis_phase = true;
analysis_shufcorrectincorrect = false; % this only matter is analysis_phase is true

test_stats = false; % this flag will trigger a downsampling procedure for the number of repetitions used in the O-score Z-scoring
plot_lastsubj = false;

% randomizations
% if testing the # of repetitions for O-score Z-scoring:
tnrep = [50:50:1000,1100:100:2000];% subsamplings
% if doing phase shuffling analyses
nshuffle_ds = 100;
nshuffle_corincor = 500;

% for stats
alpha = 0.05;
Zthres = norminv(1-alpha,0,1);

%% set configs

% time step for ACH
dt = 0.001; %s

% oscore
cfg_oscore = [];
cfg_oscore.quantlim = [0.05,0.95];
cfg_oscore.fs = 1/dt;
cfg_oscore.mincycles = 3;
cfg_oscore.smoothach = true;%
cfg_oscore.smoothwind = 0.002; %s
cfg_oscore.peakwind = 0.008; %s
cfg_oscore.taper = 'hanning';
cfg_oscore.warnings = 'off';
cfg_oscore.fcor = false;
cfg_oscore.nrep = 500;
cfg_oscore.flim = [0.5,40];
cfg_oscore.keep_trend = true;
cfg_oscore.trend_dist = {'Gamma'};
cfg_oscore.trend_ddt = 0.0005; %s for generating random data
cfg_oscore.trend_alpha = alpha;
cfg_oscore.plot = false;
if test_stats
    cfg_oscore.nrep = tnrep(end);
end

% to create time traces
cfg_mct = [];
cfg_mct.sd_smooth = [];
cfg_mct.dt = dt;
cfg_mct.removeVal = 0;
cfg_mct.warnings = 'off';

% for filtering and hilbert transform
cfg_hilb = [];
cfg_hilb.tol = 10;
cfg_hilb.fs = 1/dt;
cfg_hilb.freq_filt = '[fosc_ref-0.5, fosc_ref+0.5]';
cfg_hilb.smoothwind = '1/(fosc_ref*8)';

% some colors
cl = [];
cl.blue = [27,117,187]/255;
cl.purple = [0.4 0 0.5];
cl.green = [0.1 0.5 0.2];
cl.gray = [0.6,0.6,0.6];
cl.white = [1,1,1];

%% set relevant task data

switch taskphase
    
    case 'visual'
        experiments = 1:4;
        
    case 'encoding'
        experiments = [5:11,1213];
        
    case 'retrieval'
        experiments = [10,11,1213];
        
    case 'catch'
        experiments = [10,11,1213];
        
    case 'ret&catch'
        experiments = 5:9;
        
    otherwise
        error('Unknown task phase specified')
end

%% load data

dataBeh = table;
for e = experiments
    load([datapath, 'BehavioralData_Experiment', num2str(e),'.mat'])
    dataBeh = [dataBeh;dataBehavior];
end
clear dataBehavior

subj = unique(dataBeh.SubjID);
nsubj = length(subj);

%% pre-allocate memory

subjID = cell(nsubj,1);
subjIncl = ones(nsubj,1);
acc = zeros(nsubj,1);

% RT data for plotting
[RT_correct_st, RT_incorrect_st] = deal(cell(nsubj,1));

% oscore analyses
[OscoreZ_correct, Oscore_correct_st, fosc_correct_st] = deal(nan(nsubj,1));
OscoreZ_correct_test = nan(nsubj,length(tnrep),50);

% phase analyses
[phase_correct_st, phase_incorrect_st,...
    phase_correct_shuf_st, phase_incorrect_shuf_st] = deal(cell(nsubj,1));
[pval_cor,v_cor,pval_incor,v_incor] = deal(nan(nsubj,nshuffle_ds));
subj_downs = 0;
subj_phase = 0;

%% compute oscillation score 

for ns = 1:nsubj
    fprintf('Processing subject %i of %i: ', ns, nsubj);
    
    cfgt_oscore = cfg_oscore;
    cfgt_mct = cfg_mct;
    
    % pick data from this subject
    data = dataBeh(dataBeh.SubjID == subj(ns) ...
        & strcmp(dataBeh.TrialType,taskphase),:);
    
    subjID{ns} = subj(ns);
    
    % find performance
    Ntrials = height(data);
    Ncor = sum(data.Acc == 1);
    Nincor = sum(data.Acc == 0);
    Ntimeout = sum(isnan(data.Acc));
    
    % find reaction times
    RT_correct = data.RTs(data.Acc == 1)';
    RT_incorrect = data.RTs(data.Acc == 0)';
    
    % store for plotting
    RT_correct_st{ns} = RT_correct;
    RT_incorrect_st{ns} = RT_incorrect;
    
    acc(ns) = Ncor / (Ncor+Nincor);
    % test performance
    if Ncor < binoinv(0.95,Ncor+Nincor, 1/2)
        subjIncl(ns) = 0;
        fprintf('Subject %i excluded for bad performance\n', ns)
        continue
    elseif length(RT_correct)<10
        subjIncl(ns) = 0;
        fprintf('Subject %i excluded due to too few data points\n', ns)
        continue
    end
    
    % make continuous time lines for ACH
    [sign_correct, tspan_correct] = makeContinuousTrace(cfgt_mct,RT_correct);
    
    % oscillation score - correct trials
    cfgt_oscore.fpeak = [];
    [Oscore_correct, fosc_correct,~,~,~] = oscillationScore(cfgt_oscore, sign_correct);
    Oscore_correct_st(ns) = Oscore_correct;
    fosc_correct_st(ns) = fosc_correct;
    
    % O-score stats
    cfgt_oscore.fpeak = fosc_correct;
    [Oscore_correct_rep, fosc_correct_rep,~,trendfit] = oscillationScoreStats(cfgt_oscore, sign_correct);
    
    % compute Z-score
    OscoreZ_correct(ns) = (log(Oscore_correct)' - nanmean(log(Oscore_correct_rep))) ./ nanstd(log(Oscore_correct_rep));
    
    % if testing the # repetitions for O-score stats
    if test_stats
        for nz = 1:length(tnrep)
            if length(Oscore_correct_rep) < cfgt_oscore.nrep
                continue
            end
            for rr = 1:50
                dum = Oscore_correct_rep(randperm(cfgt_oscore.nrep, tnrep(nz)));
                OscoreZ_correct_test(ns,nz,rr) = (log(Oscore_correct)' - nanmean(log(dum))) ./ nanstd(log(dum));
            end
        end
    end
    
    if analysis_phase && OscoreZ_correct(ns)>=Zthres % only include subjects with significant Oscore
        
        subj_phase = subj_phase+1;
        
        refdumcor = RT_correct;
        dumcor = RT_correct;
        
        fosc_ref = fosc_correct_st(ns);
        
        % correct trials - the leave one out way
        for i = 1:length(dumcor)
            
            [b,id] = ismember(dumcor(i),refdumcor);
            if b
                refdumdumcor = refdumcor([1:id-1,id+1:end]);
            else
                refdumdumcor = refdumcor;
            end
            
            dumst = cfgt_mct.sd_smooth;
            cfgt_mct.sd_smooth = eval(cfg_hilb.smoothwind); %1/(fosc_ref*8);%cfgt_oscore.smoothwind_phase;
            [signsm_all, tspansm_all] = makeContinuousTrace(cfgt_mct,refdumdumcor);
            cfgt_mct.sd_smooth = dumst;
            
            cfgt_hilb = cfg_hilb;
            cfgt_hilb.freqlim = eval(cfg_hilb.freq_filt); %[fosc_ref-0.5, fosc_ref+0.5];
            signsm_all_hilb = narrowbandHilbert(cfgt_hilb,signsm_all);
            phase_all_hilb = angle(signsm_all_hilb);
            
            [~,id] = ismember(round(dumcor(i)/dt),round(tspansm_all/dt));
            if id == 0 || dumcor(i) < min(refdumdumcor)
                % we cannot compute phases for the first and last
                % responses, as there is no reliable phase trace around
                % them
                phase_correct_st{ns}(i) = NaN;
            else
                phase_correct_st{ns}(i) = phase_all_hilb(id);
            end
        end
        
        % incorrect trials
        dumincor = RT_incorrect;
        
        cfgt_mct.sd_smooth = eval(cfg_hilb.smoothwind); 
        [signsm_all, tspansm_all] = makeContinuousTrace(cfgt_mct,refdumcor);
        cfgt_mct.sd_smooth = dumst;
        
        cfgt_hilb = cfg_hilb;
        cfgt_hilb.freqlim = eval(cfg_hilb.freq_filt);
        signsm_all_hilb = narrowbandHilbert(cfgt_hilb,signsm_all);
        phase_all_hilb = angle(signsm_all_hilb);
        
        [~,ids] = ismember(round(dumincor/dt),round(tspansm_all/dt));
        ids = ids(ids>0 & dumincor>min(refdumcor));
        phase_incorrect_st{ns} = phase_all_hilb(ids);

        
        % label shuffle to test phase modulation of correct vs incorrect
        if analysis_shufcorrectincorrect
            phase_correct_shuf_dum = nan(nshuffle_corincor,length(dumcor));
            phase_incorrect_shuf_dum = nan(nshuffle_corincor,length(dumincor));
            
            fosc_ref = fosc_correct_st(ns);
            
            for r = 1:nshuffle_corincor
                
                % shuffle correct & incorrect
                RT_all = [dumcor, dumincor];
                randids = randperm(length(RT_all),length(dumincor));
                RT_incorrect_shuf = RT_all(randids);
                RT_correct_shuf = RT_all(setdiff(1:length(RT_all),randids));
                
                refdumcorshuf = RT_correct_shuf;
                
                dumst = cfgt_mct.sd_smooth;
                cfgt_mct.sd_smooth = eval(cfg_hilb.smoothwind);
                [signsm_correct_shuf, tspansm_correct_shuf] = makeContinuousTrace(cfgt_mct,refdumcorshuf);
                cfgt_mct.sd_smooth = dumst;
                
                % filter in narrow band around the peak freq and do hilbert transfort
                cfgt_hilb = cfg_hilb;
                cfgt_hilb.freqlim = eval(cfg_hilb.freq_filt);
                signsm_hilb_shuf = narrowbandHilbert(cfgt_hilb, signsm_correct_shuf);
                phase_hilb_shuf = angle(signsm_hilb_shuf);
                
                [~,ids] = ismember(round(RT_incorrect_shuf/dt),round(tspansm_correct_shuf/dt));
                ids = ids(ids>0 & RT_incorrect_shuf > min(RT_correct_shuf));
                phase_incorrect_shuf_dum(r,1:length(ids)) = phase_hilb_shuf(ids);
                
                for i = 1:length(RT_correct_shuf)
                    
                    [b,id] = ismember(RT_correct_shuf(i),refdumcorshuf);
                    if b
                        refdumdumcor = refdumcorshuf([1:id-1,id+1:end]);
                    else
                        refdumdumcor = refdumcorshuf;
                    end
                    
                    dumst = cfgt_mct.sd_smooth;
                    cfgt_mct.sd_smooth = eval(cfg_hilb.smoothwind);
                    [signsm_correct_shuf, tspansm_correct_shuf] = makeContinuousTrace(cfgt_mct,refdumdumcor);
                    cfgt_mct.sd_smooth = dumst;
                    
                    % filter in narrow band around the peak freq and do hilbert transfort
                    cfgt_hilb = cfg_hilb;
                    cfgt_hilb.freqlim = eval(cfg_hilb.freq_filt);
                    signsm_hilb_shuf = narrowbandHilbert(cfgt_hilb, signsm_correct_shuf);
                    phase_hilb_shuf = angle(signsm_hilb_shuf);
                    
                    [~,id] = ismember(round(RT_correct_shuf(i)/dt),round(tspansm_correct_shuf/dt));
                    if id == 0 || RT_correct_shuf(i)<min(refdumdumcor)
                        phase_correct_shuf_dum(r,i) = NaN;
                    else
                        phase_correct_shuf_dum(r,i) = phase_hilb_shuf(id);
                    end
                end
            end
            phase_incorrect_shuf_st{ns} = phase_incorrect_shuf_dum;
            phase_correct_shuf_st{ns} = phase_correct_shuf_dum;
        end
        
        % subsample correct trials per subject
        if length(phase_incorrect_st{ns}) >= 10
            
            fosc_ref = fosc_correct_st(ns);
            
            % count the number of included subjects
            subj_downs = subj_downs + 1;
            
            % downsample per subject
            if length(dumincor) < length(refdumcor)
                for i = 1:nshuffle_ds
                    subids = randperm(length(refdumcor),length(dumincor));
                    othids = setdiff(1:length(refdumcor),subids);
                    
                    dumst = cfgt_mct.sd_smooth;
                    cfgt_mct.sd_smooth = eval(cfg_hilb.smoothwind);
                    [signsm_correct_sub, tspansm_correct_sub] = makeContinuousTrace(cfgt_mct,refdumcor(othids));
                    cfgt_mct.sd_smooth = dumst;
                    
                    cfgt_hilb = cfg_hilb;
                    cfgt_hilb.freqlim = eval(cfg_hilb.freq_filt);
                    signsm_hilb_sub = narrowbandHilbert(cfgt_hilb, signsm_correct_sub);
                    phase_hilb_sub = angle(signsm_hilb_sub);
                    
                    [~,ids] = ismember(round(refdumcor(subids)/dt),round(tspansm_correct_sub/dt));
                    ids = ids(ids>0);
                    phase_correct_sub = phase_hilb_sub(ids);
                    
                    [~,ids] = ismember(round(dumincor/dt),round(tspansm_correct_sub/dt));
                    ids = ids(ids>0);
                    phase_incorrect_sub = phase_hilb_sub(ids);
                    
                    [pval_cor(ns,i),v_cor(ns,i)] = circ_vtest(phase_correct_sub(~isnan(phase_correct_sub)),0);
                    [dum1,dum2] = circ_vtest(phase_incorrect_sub(~isnan(phase_incorrect_sub)),0);
                    pval_incor(ns,:) = dum1;
                    v_incor(ns,:) = dum2;
                    
                end
            else
                % this situation shouldn't occur, given the earlier test
                % for performance, but just in case: exclude this subject.
                
                pval_cor(ns,:) = NaN;
                v_cor(ns,:) = NaN;
                pval_incor(ns,:) = NaN;
                v_incor(ns,:) = NaN;
            end
        end
    end
    
    fprintf('done!\n')
end

fprintf('All subjects completed!\n')

%% save data

if save_data
    save([storepath, 'OscoreAnalysis_',taskphase, storename,'.mat'], ...
        'cfg_oscore','cfg_mct',...
        'subjID', 'subjIncl', 'acc',...
        'RT_correct_st', 'RT_incorrect_st',...
        'OscoreZ_correct', 'Oscore_correct_st', 'fosc_correct_st')
    if analysis_shufcorrectincorrect
        save([storepath, 'PhaseAnalysis_',taskphase, storename,'.mat'], ...
            'cfg_hilb', 'subj_downs','subj_phase',...
            'phase_correct_st', 'phase_incorrect_st',...
            'pval_cor','v_cor','pval_incor','v_incor',...
            'phase_correct_shuf_st','phase_incorrect_shuf_st')
    end
    if test_stats
        save([storepath, 'NrepAnalysis_',taskphase, storename,'.mat'], ...
            'OscoreZ_correct_test')
    end
end

%% plot test stats (Figure S2a)

if test_stats
    
    figure;
    subplot(211); hold on;
    dat = repmat(OscoreZ_correct,[1,length(tnrep),50]);
    errors = reshape(permute(...
        abs(OscoreZ_correct_test - dat) ./ dat ...
        ,[2,1,3]),...
        [length(tnrep),nsubj*50]);
    id = ~isnan(errors);
    plot([0,2000],[5,5],'k:')
    xdat = repmat(tnrep',[1,size(errors,2)]);
    simplebox(xdat(id), 100*errors(id), [0,0,0]);
    ylim([0,40])
    xlabel('# permutations')
    ylabel('Error (%)')
    title('O-score (Z-scored)')
    
    subplot(212); hold on;
    dum = squeeze(nansum(abs((OscoreZ_correct_test >= Zthres) - ...
        repmat(OscoreZ_correct >= Zthres,[1,length(tnrep),50])),1) / (nsubj));
    plot([0,2000],[5,5],'k:')
    xdat = repmat(tnrep',[1,size(dum,2)]);
    simplebox(xdat(:), 100*dum(:),[0,0,0]);
    ylim([0,20])
    xlabel('# permutations')
    ylabel('Error (%)')
    title('Detecting significance')
    
    if save_figures
        saveas(gcf,[storepath, 'TestNrep_',taskphase, '_',catchtype,storename,'.fig'])
    end
    
end

%% plot last subject (Figure 2a)
% plots the subject run last as an example

if plot_lastsubj
    
    cfg = cfg_mct;
    cfg.sd_smooth = 0.050;
    [signsm_correct, tspansm_correct] = makeContinuousTrace(cfg,RT_correct);
    
    figure;
    subplot(3,1,1:2);hold on;
    
    % smoothed
    plot(cfg_mct.dt:cfg_mct.dt:(length(trendfit.trace)*cfg_mct.dt),trendfit.trace*length(RT_correct), 'color',cl.gray, 'linewidth',2)
    plot(tspansm_correct,signsm_correct, 'color',cl.green, 'linewidth',2)
    
    % scatter
    dum = max(signsm_correct);
    p1 = scatter(RT_correct,dum+0.025+0.015*rand(size(RT_correct)),15, cl.green, 'filled');
    p2 = scatter(RT_incorrect,dum+0.005+0.005*rand(size(RT_incorrect)),15, cl.purple, 'filled');
    title(['Subject # ', num2str(ns), '; O-scoreZ = ', num2str(OscoreZ_correct(ns),3), '; fosc = ', num2str(fosc_correct_st(ns))])
    legend([p1,p2], {'Correct','Incorrect'})
    xlabel('Time (s)')
    ylabel('Response density (a.u.)')
    xlim([0,8])
    
    subplot(313);hold on;
    plot([0,8], [0,0],'k')
    plot(cfg_mct.dt:cfg_mct.dt:(length(trendfit.trace)*cfg_mct.dt),signsm_correct - trendfit.trace*length(RT_correct), 'color',cl.green, 'linewidth',2)
    xlabel('Time (s)')
    ylabel('Rel. resp. density (a.u.)')
    xlim([0,8])
    
    if save_figures
        saveas(gcf,[storepath, 'ExampleSubject_',taskphase, '_',catchtype,storename,'.fig'])
    end
end

%% plot oscillation score (Figure 3 and S3)

% plot Oscore for correct trials
fmax = 30;
figure('Position',[50,50,800,600]);
subplot(2,4,[1,5]); hold on
dum = zeros(size(OscoreZ_correct));
dum(OscoreZ_correct>=Zthres) = 1;
fcolors = [cl.gray;cl.white;cl.blue];
lcolors = [cl.gray;cl.blue;cl.blue];
dat = [sum(subjIncl==0),sum(dum==0)-sum(subjIncl==0),sum(dum==1)];
for i = 1:3
    bar(i, dat(i), 'facecolor', fcolors(i,:), 'edgecolor',lcolors(i,:))
end
set(gca,'xtick',1:3,'xticklabel',{'Excl.','No osc.','Osc'})
rotateXLabels(gca,90)
xlim([0,4])
ylabel('# subjects')

% second level stats
dum = OscoreZ_correct - Zthres;
vargroup = nansum((dum - nanmean(dum)).^2)/(sum(subjIncl)-1);
Tpop_Os = (nanmean(dum)) / sqrt(vargroup/(sum(subjIncl)-1));
title({['T = ', num2str(Tpop_Os), '; p = ',num2str(1-tcdf(Tpop_Os,sum(subjIncl)-1),2)]})
nanmean(OscoreZ_correct)

subplot(2,4,2:4); hold on;
scatter(fosc_correct_st(OscoreZ_correct<Zthres),OscoreZ_correct(OscoreZ_correct<Zthres), 20, cl.blue)
scatter(fosc_correct_st(OscoreZ_correct>=Zthres),OscoreZ_correct(OscoreZ_correct>=Zthres), 20, cl.blue, 'filled')
plot([0,fmax],[Zthres,Zthres],'k:')
plot([0,fmax],[-Zthres,-Zthres],'k:')
xlim([0,fmax])
ylim([-2,15])
title('O-score (Z-scored)')

subplot(2,4,6:8)
hi = hist(fosc_correct_st(OscoreZ_correct<Zthres),0:fmax);
hc = hist(fosc_correct_st(OscoreZ_correct>=Zthres),0:fmax);
b = bar(0:fmax,[hi;hc]', 'edgecolor',cl.blue);
b(1).FaceColor = cl.white;
b(2).FaceColor = cl.blue;
xlabel('Peak frequency (Hz)')
title('Frequency')
xlim([0,fmax])

if save_figures
    saveas(gcf,[storepath, 'Oscore_',taskphase, '_',catchtype,storename,'.fig'])
end

% scatter plot of Oscores (Z-scored and raw)
figure;
subplot(121); hold on
plot([0,2],[0,0],'k')
plot([0,2],[Zthres,Zthres],'k:')
scatter(ones(size(OscoreZ_correct))-0.2*rand(size(OscoreZ_correct)), OscoreZ_correct,10,cl.blue)
simplebox(0.2+ones(size(OscoreZ_correct)), OscoreZ_correct, cl.blue)
set(gca,'xtick',1,'xticklabel',taskphase)
ylabel('O-score (T-scored)')
ylim([-1,12])

subplot(122); hold on
plot([0,2],[0,0],'k')
scatter(ones(size(Oscore_correct_st))-0.2*rand(size(Oscore_correct_st)), Oscore_correct_st,10,cl.blue)
simplebox(0.2+ones(size(Oscore_correct_st)), Oscore_correct_st, cl.blue)
ylabel('O-score')
set(gca,'xtick',1,'xticklabel',taskphase)
ylim([-10,80])

if save_figures
    saveas(gcf,[storepath, 'Oscore_scatter_',taskphase, '_',catchtype,storename,'.fig'])
end

%% plot phase analyses (Figure 4)

if analysis_phase
    
    % plot V-stats for subsampling
    figure; subplot(121); hold on;
    dum1 = 0.2*rand(size(v_cor,1),1)-0.1;
    dum2 = 0.2*rand(size(v_incor,1),1)-0.1;
    
    % - plot means per subject
    plot([dum1+1,dum2+2]',[nanmean(v_cor,2),nanmean(v_incor,2)]','k')%
    scatter(dum1+ones(size(v_cor,1),1), nanmean(v_cor,2),20,cl.green, 'filled')
    scatter(dum2+2*ones(size(v_incor,1),1), nanmean(v_incor,2), 20,cl.purple, 'filled')
    
    % - plot entire distibution
    simplebox(0.7*ones(size(v_cor)), v_cor,cl.green)
    simplebox(2.3*ones(size(v_incor)), v_incor,cl.purple)
    ylim([-15,15])
    xlim([0.5,2.5])
    ylabel('V-stat')
    set(gca,'xtick' ,[1,2],'xticklabel' ,{'Correct','Incorrect'})
    
    [h,p1,~, stats] = ttest(nanmean(v_cor,2),nanmean(v_incor,2));
    title(['p = ' ,num2str(p1,3)])
    
    for i = 1:nshuffle_ds
        [h(i),p(i)] = ttest(v_cor(:,i),v_incor(:,i));
    end
    pval = 1-nanmean(h);
    suptitle(['Perm. stats: p-val = ', num2str(pval), '; Nsubj = ', num2str(subj_downs)])
    CI = quantile(p,[0.025,0.975]);
    title(['p = ' ,num2str(p1,3), '; CI = [',num2str(CI),']'])
    
    subplot(122); hold on;
    b = bar(1,sum(pval_cor(:)<0.05)/sum(~isnan(pval_cor(:))));
    b(1).FaceColor = cl.green;
    b = bar(2,sum(pval_incor(:)<0.05)/sum(~isnan(pval_incor(:))));
    b(1).FaceColor = cl.purple;
    set(gca,'xtick' ,[1,2],'xticklabel' ,{'Correct','Incorrect'})
    ylabel('# subjects with significant V-stat')
    ylim([0,0.5])
    
    if save_figures
        saveas(gcf,[storepath, 'PhaseModulation_CorrectIncorrect_',taskphase, '_',catchtype,storename,'.fig'])
    end
    
    % plot phase distributions across subjects
    nbins = 10; % number of bins for phase plotting
    
    figure; hold on
    ls = linspace(-pi+pi/nbins,pi-pi/nbins,nbins);
    dumc = cat(2,phase_correct_st{:});
    dumi = cat(1,phase_incorrect_st{:});
    hc = hist(dumc(~isnan(dumc)),ls)/ (length(dumc(~isnan(dumc)))/nbins);
    hi = hist(dumi(~isnan(dumi)),ls) / (length(dumi(~isnan(dumi)))/nbins);
    plot([0,0],[-50,50],'k: ')
    plot([-pi,pi],[0,0],'k')
    plot(ls,hi*100-100, 'o-', 'color',cl.purple, 'linewidth',3, 'markerfacecolor',cl.purple)
    plot(ls,hc*100-100, 'o-','color', cl.green, 'linewidth',3, 'markerfacecolor',cl.green)
    set(gca,'xtick', [-pi,-pi/2,0,pi/2,pi], 'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
    ylabel('% change relative to uniform')
    xlabel('Phase of response')
    ylim([-50,50])
    xlim([-pi,pi])
    text(0.9,27,['Correct trials; p = ', num2str(3*circ_vtest(dumc(~isnan(dumc)),0),3)], 'color',cl.green)
    text(0.9,24,['Incorrect trials; p = ', num2str(3*circ_vtest(dumi(~isnan(dumi)),0),3)],'color',cl.purple)
    title('Original data')
        
    if save_figures
        saveas(gcf,[storepath, 'Phase_CorrectIncorrectLinAx_',taskphase, '_',catchtype,storename,'.fig'])
    end
    
    if analysis_shufcorrectincorrect
        % correct & incorrect phase distributions - with permutation stats
        
        % shuffled data alone
        figure; hold on
        ls = linspace(-pi+pi/nbins,pi-pi/nbins,nbins);
        datshi = cat(2,phase_incorrect_shuf_st{:});
        datshi = datshi(~isnan(datshi));
        hi = hist(datshi(:),ls)/ (length(datshi(:))/nbins);
        datshc = cat(2,phase_correct_shuf_st{:});
        datshc = datshc(~isnan(datshc));
        hc = hist(datshc(:),ls) / (length(datshc(:))/nbins);
        plot([0,0],[-50,50],'k: ')
        plot([-pi,pi],[0,0],'k')
        plot(ls,hi*100-100, 'o-', 'color',cl.purple, 'linewidth',3, 'markerfacecolor',cl.purple)
        plot(ls,hc*100-100, 'o-','color', cl.green, 'linewidth',3, 'markerfacecolor',cl.green)
        set(gca,'xtick', [-pi,-pi/2,0,pi/2,pi], 'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
        ylabel('Fraction of responses')
        xlabel('Phase of response')
        ylim([-30,30])
        xlim([-pi,pi])
        text(0.9,27,['Correct trials; p = ', num2str(3*circ_vtest(datshc(~isnan(datshc(:))),0),3)], 'color',cl.green)
        text(0.9,24,['Incorrect trials; p = ', num2str(3*circ_vtest(datshi(~isnan(datshi(:))),0),3)],'color',cl.purple)
        title('Shuffled data')
        
        if save_figures
            saveas(gcf,[storepath, 'Phase_CorrectIncorrectLinAx_shuffled_',taskphase, '_',catchtype,storename,'.fig'])
        end
        
        % compute permutation stats
        % - observed data
        [p,v_co] = circ_vtest(dumc(~isnan(dumc)),0.0);
        [p,v_io] = circ_vtest(dumi(~isnan(dumi)),0.0);
        amp_orig = v_co - v_io;
        
        % - shuffled data
        datcst = cat(2,phase_correct_shuf_st{:});
        datist = cat(2,phase_incorrect_shuf_st{:});
        amp_shuf = zeros(nshuffle_corincor,1);
        for r = 1:nshuffle_corincor
            dum = datcst(r,:);
            [~,v_c(r)] = circ_vtest(dum(~isnan(dum)),0);
            dum = datist(r,:);
            [~,v_i(r)] = circ_vtest(dum(~isnan(dum)),0);
            amp_shuf(r) = v_c(r)- v_i(r);
        end
        
        figure; hold on
        ls = linspace(-pi+pi/nbins,pi-pi/nbins,nbins);
        dumc = cat(2,phase_correct_st{:});
        dumi = cat(1,phase_incorrect_st{:});
        hc = hist(dumc(~isnan(dumc)),ls)/ (length(dumc(~isnan(dumc)))/nbins);
        hi = hist(dumi(~isnan(dumi)),ls) / (length(dumi(~isnan(dumi)))/nbins);
        plot([0,0],[-50,50],'k: ')
        plot([-pi,pi],[0,0],'k')
        plot(ls,hi*100-100, 'o-', 'color',cl.purple, 'linewidth',3, 'markerfacecolor',cl.purple)
        plot(ls,hc*100-100, 'o-','color', cl.green, 'linewidth',3, 'markerfacecolor',cl.green)
        set(gca,'xtick', [-pi,-pi/2,0,pi/2,pi], 'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
        ylabel('% change relative to uniform')
        xlabel('Phase of response')
        ylim([-50,50])
        xlim([-pi,pi])
        text(0.9,27,['Correct trials; p = ', num2str(3*circ_vtest(dumc(~isnan(dumc)),0),3)], 'color',cl.green)
        text(0.9,24,['Incorrect trials; p = ', num2str(3*circ_vtest(dumi(~isnan(dumi)),0),3)],'color',cl.purple)
        title('Original data')
        
        % test difference of V-stats:
        % H1: original dist correct-incorrect is more phase modulated than shuffled data (V_orig>V_shuf)
        % test whether orig is significantly LESS modulated than shuffled
        pval = sum(amp_shuf>amp_orig)/length(amp_shuf);
        
        title(['V-stat permutation test: p = ', num2str(pval), ' Nsubj = ', num2str(subj_phase)])
        
        if save_figures
            saveas(gcf,[storepath, 'Phase_CorrectIncorrectLinAx_PermStats_',taskphase, '_',catchtype,storename,'.fig'])
        end
    end
end
