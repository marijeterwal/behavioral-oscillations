
%{
This script analyses the correct versus incorrect trials contrast on PPC data
and produces plots for Figures 5 and S4 from:

Ter Wal, M., Linde Domingo, J., Lifanov, J., Roux, F., Kolibius, L., 
Gollwitzer, S., Lang, J., Hamer, H., Rollings, D., Sawlani, V., 
Chelvarajah, R., Staresina, B., Hanslmayr, S., Wimber, M. 
Theta rhythmicity governs the timing of behavioural and hippocampal 
responses in humans specifically during memory- dependent tasks. 
bioRxiv 2020.11.09.374264; DOI: https://doi.org/10.1101/2020.11.09.374264

The data are available here:
Ter Wal, M., Linde Domingo, J., Lifanov, J., Roux, F., Kolibius, L., 
Gollwitzer, S., Lang, J., Hamer, H., Rollings, D., Sawlani, V., Chelvarajah, 
R., Staresina, B., Hanslmayr, S., Wimber, M., (2020), Data for: Theta 
rhythmicity governs the timing of behavioural and hippocampal responses in 
humans specifically during memory- dependent tasks. figshare. Collection. 
DOI: 10.6084/m9.figshare.c.5192567


First, the script loads all PPC data and computes the contrast between correct
and incorrect trials. It performs cluster detection on both the original
data and the trial-shuffled PPC data.

% The cluster-based analysis consists of the following steps:
% step 1: T-score correct vs incorrect trials for the original data
% step 2: find 2D clusters in original data
% step 3: T-score shuffled-correct vs shuffled-incorrect trials for each repetition
% step 4: find 2D clusters in shuffled data
% step 5: compare cluster T-stats

Finally, the code plots the T-contrasts of the original data and
statistics, as well as the reaction time distributions.

Note: this script requires loading of a lot of data, which makes it slow
and memory intensive!

% The script relies on the following functions:
- UScoreMat
- pairedTscoreMat
- detectClusters and dependencies
- nonParamPVal

External:
- brewermap

Marije ter Wal - 2020
m.j.terwalm@bham.ac.uk

%}

clear all

%% settings

datapath = 'YourPathGoesHere';
datapathRT = 'YourPathGoesHere';

taskphases = {'encoding', 'retrieval','catch'};
nevents = [3,2,2];

trialtypes = {'correct','incorrect'};
acc = [1,0];

nrep = 100;

alpha = 0.05; % for cluster detection threshold
alphac = 0.05; % for permutation test

% colors for plotting
col = [];
col.red = [190 30 45]/255;
col.blue = [27,117,187]/255;
col.purple = [0.4 0 0.5];
col.green = [0.1 0.6 0.2];
col.gray = [0.6,0.6,0.6];

% cluster detection
clus_method = 'sumZscore';

cfg = [];
cfg.method = 'extract';
cfg.freqweight = 1.5;
cfg.timeweight = 1;

cfg.splitClusters = true;
cfg.peakHeight = 0.05:0.05:0.9;
cfg.minNewClusSize = 0.05;

%% find subject IDs

files = dir([datapath,'PPC_SubjID-1213*_encoding_correct.mat']);
subjIDs = arrayfun(@(x) str2num(files(x).name(12:17)),1:length(files));
nsubj = length(subjIDs);

assert(nsubj>0, 'No data found')

%% load original PPC data

% preallocate
for tp = 1:length(taskphases)
    for tr = 1:length(trialtypes)
        eval(['PPCZ_',taskphases{tp},'_',trialtypes{tr},' = cell(nsubj,nevents(tp));']);
    end
end
separateCue = ones(nsubj,1);

for s = 1:nsubj
    
    subjID = num2str(subjIDs(s));
    fprintf('Loading original data of subject %i of %i\n', s, nsubj)
    
    for tp = 1:length(taskphases)
        
        for tr = 1:length(trialtypes)
            % encoding
            load([datapath, 'PPC_SubjID-',subjID,'_',taskphases{tp},'_',trialtypes{tr},'.mat'])
            
            for ev = 1:length(PPC.PPC)
                % normalize every channel to its baseline
                PPCZ = nan(size(PPC.PPC{ev}));
                for ch = 1:size(PPC.PPC{ev},1)
                    ref = repmat(permute(squeeze(PPC.PPC_baseline{1}(ch,:,:)),[2,1]),[1,1,size(PPC.PPC{ev},3)]);
                    [PPCZ(ch,:,:),~] = UScoreMat(squeeze(PPC.PPC{ev}(ch,:,:)), ref, 0.05);
                end
                
                % average within bundles
                for b = unique(PPC.bundle)'
                    bids = find(PPC.bundle==b);
                    eval(['PPCZ_',taskphases{tp},'_',trialtypes{tr},'{s,ev}(b,:,:) = squeeze(nanmean(PPCZ(bids,:,:),1));'])
                end
            end
            eval(['time_',taskphases{tp},' = PPC.time;']);
            freq = PPC.freq;
            dt = 1/PPC.fsample;
            
            clear PPC PPCZ ref
        end
        
        if strcmp(taskphases{tp},'encoding') && all(PPCZ_encoding_correct{s,1}(1,1,:)==PPCZ_encoding_correct{s,2}(1,1,:),3)
            separateCue(s) = 0; 
        end
    end
end

subjectSelCue = find(separateCue);
subjectSel = 1:nsubj;

%% load label-shuffled PPC data

% preallocate
for tp = 1:length(taskphases)
    for tr = 1:length(trialtypes)
        eval(['PPCZ_',taskphases{tp},'_',trialtypes{tr},'_shuf = cell(nsubj,nrep,nevents(tp));']);
    end
end

for s = 1:nsubj
    
    subjID = num2str(subjIDs(s));
    fprintf('Loading shuffled data of subject %i of %i\n', s, nsubj)
    
    for tp = 1:length(taskphases)
        
        for tr = 1:length(trialtypes)
            % encoding
            load([datapath, 'PPC_SubjID-',subjID,'_',taskphases{tp},'_',trialtypes{tr},'_trialshuffle.mat'])
            
            for nr = 1:nrep
                for ev = 1:length(PPC.PPC)
                    % normalize every channel to its baseline
                    sz = size(PPC.PPC{ev});
                    PPCZ = nan(sz(1:3));
                    for ch = 1:size(PPC.PPC{ev},1)
                        ref = repmat(permute(squeeze(PPC.PPC_baseline{1}(ch,:,:,nr)),[2,1]),[1,1,size(PPC.PPC{ev},3)]);
                        [PPCZ(ch,:,:),~] = UScoreMat(squeeze(PPC.PPC{ev}(ch,:,:,nr)), ref, 0.05);
                    end
                    
                    % average within bundles
                    for b = unique(PPC.bundle)'
                        bids = find(PPC.bundle==b);
                        eval(['PPCZ_',taskphases{tp},'_',trialtypes{tr},'_shuf{s,nr,ev}(b,:,:) = squeeze(nanmean(PPCZ(bids,:,:),1));'])
                    end
                end
            end
            eval(['time_',taskphases{tp},' = PPC.time;']);
            freq = PPC.freq;
            
            clear PPC PPCZ ref
        end
    end
end

%% STEP 1: compute second level T-score for original data

for tp = 1:length(taskphases)
    for ev = 1:nevents(tp)
        if strcmp(taskphases{tp},'encoding') && ev == 1
            dat1 = eval(['cat(1,PPCZ_',taskphases{tp},'_',trialtypes{1},'{subjectSelCue,ev})']);
            dat2 = eval(['cat(1,PPCZ_',taskphases{tp},'_',trialtypes{2},'{subjectSelCue,ev})']);
        else
            dat1 = eval(['cat(1,PPCZ_',taskphases{tp},'_',trialtypes{1},'{subjectSel,ev})']);
            dat2 = eval(['cat(1,PPCZ_',taskphases{tp},'_',trialtypes{2},'{subjectSel,ev})']);
        end
       eval(['[t_',taskphases{tp},num2str(ev),',nu_',taskphases{tp},num2str(ev),'] = pairedTscoreMat(dat1,dat2,1);'])
    end
end

%% STEP 2: find 2D clusters

% for cluster detection:
cfg.freqs = freq;
cfg.cw = 4;
cfg.dt = dt;
cfg.fs = 1/cfg.dt;

for tp = 1:length(taskphases)
    for ev = 1:nevents(tp)
        % define boundaries for cluster detection
        cfg.preEvent = eval(['time_', taskphases{tp},'{ev}(1)']);
        cfg.postEvent = eval(['time_', taskphases{tp},'{ev}(end)']);
        
        % find clusters!
        eval(['clusStruct_',taskphases{tp},num2str(ev),...
            ' = detectClusters(cfg, squeeze(t_',taskphases{tp},num2str(ev),...
            '),tinv(1-alpha,squeeze(nu_', taskphases{tp},num2str(ev),')),[],[]);']);
    end
end

%% STEP 3 and 4: compute T-max dist from shuffled data

% preallocate
for tp = 1:length(taskphases)
    for ev = 1:nevents(tp)
        eval(['clusStruct_r',taskphases{tp},num2str(ev),'_Pos = [];'])
        eval(['clusStruct_r',taskphases{tp},num2str(ev),'_Neg = [];'])
    end
end

for nr = 1:nrep
    fprintf('Analysing repetition %i of %i\n', nr,nrep)
    
    for tp = 1:length(taskphases)
        for ev = 1:nevents(tp)
            if strcmp(taskphases{tp},'encoding') && ev == 1
                dat1 = eval(['cat(1,PPCZ_',taskphases{tp},'_',trialtypes{1},'_shuf{subjectSelCue,nr,ev})']);
                dat2 = eval(['cat(1,PPCZ_',taskphases{tp},'_',trialtypes{2},'_shuf{subjectSelCue,nr,ev})']);
            else
                dat1 = eval(['cat(1,PPCZ_',taskphases{tp},'_',trialtypes{1},'_shuf{subjectSel,nr,ev})']);
                dat2 = eval(['cat(1,PPCZ_',taskphases{tp},'_',trialtypes{2},'_shuf{subjectSel,nr,ev})']);
            end
            eval(['[t_r',taskphases{tp},num2str(ev),',nu_r',taskphases{tp},num2str(ev),'] = pairedTscoreMat(dat1,dat2,1);'])
            
            % define boundaries for cluster detection
            cfg.preEvent = eval(['time_', taskphases{tp},'{ev}(1)']);
            cfg.postEvent = eval(['time_', taskphases{tp},'{ev}(end)']);
            
            % find clusters for this repetiiton
            eval(['tmpStruct_',taskphases{tp},num2str(ev),...
                ' = detectClusters(cfg, squeeze(t_r',taskphases{tp},num2str(ev),...
                '),tinv(1-alpha,squeeze(nu_r', taskphases{tp},num2str(ev),')),[],[]);']);
            
            % store positive and negative clusters
            eval(['clusStruct_r',taskphases{tp},num2str(ev),'_Pos',...
                '= [clusStruct_r',taskphases{tp},num2str(ev),'_Pos,',...
                'max([tmpStruct_', taskphases{tp},num2str(ev),'(:).',clus_method,'])];']);
            eval(['clusStruct_r',taskphases{tp},num2str(ev),'_Neg'...
                '= [clusStruct_r',taskphases{tp},num2str(ev),'_Neg,',...
                'min([tmpStruct_', taskphases{tp},num2str(ev),'(:).',clus_method,'])];']);
            
        end
    end    
end

for tp = 1:length(taskphases)
    for ev = 1:nevents(tp)
        eval(['clusStruct_r',taskphases{tp},num2str(ev),...
            '= [clusStruct_r',taskphases{tp},num2str(ev),'_Pos',...
            ', clusStruct_r',taskphases{tp},num2str(ev),'_Neg];'])
    end
end

%% STEP 5: non-parametric test

for tp = 1:length(taskphases)
    for ev = 1:nevents(tp)
        % compute p-value
        [h,p] = nonParamPVal(eval(['[clusStruct_',taskphases{tp},num2str(ev),...
            '(:).',clus_method,']']),eval(['clusStruct_r',taskphases{tp},num2str(ev)]),alphac);
        % store in struct
        dum = num2cell(p);
        eval(['[clusStruct_',taskphases{tp},num2str(ev),'.clus_pval] = dum{:};'])
        dum = num2cell(h);
        eval(['[clusStruct_',taskphases{tp},num2str(ev),'.clus_h] = dum{:};'])
    end
end

%% load RTs

% load RT data
load([datapathRT,'BehavioralData_Experiment1213.mat'])

for tp = 1:length(taskphases)
    for tr = 1:length(trialtypes)
        RTs = dataBehavior.RTs(strcmp(dataBehavior.TrialType,taskphases{tp}) & ...
            dataBehavior.Acc == acc(tr));
        eval(['RT_',taskphases{tp},'_',trialtypes{tr}, ...
            ' = RTs;']);
    end
end

% The cue period chosen for the PPC does not overlap with the stimulus
% period and hence has no responses:
RT_encoding_correct_cue = [];
RT_encoding_incorrect_cue = [];

%% Figure 5 and S4: plot PPC contrast with clusters and RT distributions

clear diff

cax = [-5,5];
yax = [1.0,12];
rax = [0,0.7];
map = flipud(brewermap(100,'RdBu'));

figure('Position',[50,50,1200,500]);
subplot(3,15,[1:2,16:17]); hold on;
pcolor(time_encoding{1}, freq, squeeze(t_encoding1)); shading interp
plot([0,0],[freq(1)-0.25,freq(end)+0.25],'color',col.gray); axis tight
for h = find([clusStruct_encoding1.clus_h]==1)
    dum = zeros(length(time_encoding{1}),length(freq));
    [~,iat] = ismember(clusStruct_encoding1(h).times,time_encoding{1});
    [~,iaf] = ismember(clusStruct_encoding1(h).freqs,freq);
    dum(sub2ind(size(dum), iat,iaf)) = 1;
    if clusStruct_encoding1(h).peakZscore > 0
        cluscol = map(end,:);
    else
        cluscol = map(1,:);
    end
    contour(time_encoding{1}, freq+[0.05,diff(freq)/2], dum',1, 'linewidth', 2, 'color', cluscol)
end
colormap(map);
axis xy; caxis(cax); ylim(yax);
ylabel('Frequency (Hz)')
title('Enc - cue')
subplot(3,15,[3:4,18:19]); hold on;
pcolor(time_encoding{2}, freq, squeeze(t_encoding2)); shading interp
plot([0,0],[freq(1)-0.25,freq(end)+0.25],'color',col.gray); axis tight
for h = find([clusStruct_encoding2.clus_h]==1)
    dum = zeros(length(time_encoding{2}),length(freq));
    [~,iat] = ismember(clusStruct_encoding2(h).times,time_encoding{2});
    [~,iaf] = ismember(clusStruct_encoding2(h).freqs,freq);
    dum(sub2ind(size(dum), iat,iaf)) = 1;
    if clusStruct_encoding2(h).peakZscore > 0
        cluscol = map(end,:);
    else
        cluscol = map(1,:);
    end
    contour(time_encoding{2}, freq+[0.05,diff(freq)/2], dum',1, 'linewidth', 2, 'color', cluscol)
end
axis xy; caxis(cax); ylim(yax);
set(gca,'ytick','')
title('Enc - stimulus')
subplot(3,15,[5:6,20:21]); hold on;
pcolor(time_encoding{3}, freq, squeeze(t_encoding3)); shading interp
plot([0,0],[freq(1)-0.25,freq(end)+0.25],'color',col.gray); axis tight
for h = find([clusStruct_encoding3.clus_h]==1)
    dum = zeros(length(time_encoding{3}),length(freq));
    [~,iat] = ismember(clusStruct_encoding3(h).times,time_encoding{3});
    [~,iaf] = ismember(clusStruct_encoding3(h).freqs,freq);
    dum(sub2ind(size(dum), iat,iaf)) = 1;
    if clusStruct_encoding3(h).peakZscore > 0
        cluscol = map(end,:);
    else
        cluscol = map(1,:);
    end
    contour(time_encoding{3}, freq+[0.05,diff(freq)/2], dum',1, 'linewidth', 2, 'color', cluscol)
end
axis xy; caxis(cax); ylim(yax);
set(gca,'ytick','')
xlabel('Time (s)')
title('Enc - response')

subplot(3,15,[7:8,22:23]); hold on;
pcolor(time_retrieval{1}, freq, squeeze(t_retrieval1)); shading interp
plot([0,0],[freq(1)-0.25,freq(end)+0.25],'color',col.gray); axis tight
for h = find([clusStruct_retrieval1.clus_h]==1)
    dum = zeros(length(time_retrieval{1}),length(freq));
    [~,iat] = ismember(clusStruct_retrieval1(h).times,time_retrieval{1});
    [~,iaf] = ismember(clusStruct_retrieval1(h).freqs,freq);
    dum(sub2ind(size(dum), iat,iaf)) = 1;
    if clusStruct_retrieval1(h).peakZscore > 0
        cluscol = map(end,:);
    else
        cluscol = map(1,:);
    end
    contour(time_retrieval{1}, freq+[0.05,diff(freq)/2], dum',1, 'linewidth', 2, 'color', cluscol)
end
axis xy; caxis(cax); ylim(yax);
set(gca,'ytick','')
xlabel('Time (s)')
title('Ret - cue')
subplot(3,15,[9:10,24:25]); hold on;
pcolor(time_retrieval{2}, freq, squeeze(t_retrieval2)); shading interp
plot([0,0],[freq(1)-0.25,freq(end)+0.25],'color',col.gray); axis tight
for h = find([clusStruct_retrieval2.clus_h]==1)
    dum = zeros(length(time_retrieval{2}),length(freq));
    [~,iat] = ismember(clusStruct_retrieval2(h).times,time_retrieval{2});
    [~,iaf] = ismember(clusStruct_retrieval2(h).freqs,freq);
    dum(sub2ind(size(dum), iat,iaf)) = 1;
    if clusStruct_retrieval2(h).peakZscore > 0
        cluscol = map(end,:);
    else
        cluscol = map(1,:);
    end
    contour(time_retrieval{2}, freq+[0.05,diff(freq)/2], dum',1, 'linewidth', 2, 'color', cluscol)
end
axis xy; caxis(cax); ylim(yax);
set(gca,'ytick','')
xlabel('Time (s)')
title('Ret - response')

subplot(3,15,[11:12,26:27]); hold on;
pcolor(time_catch{1}, freq, squeeze(t_catch1)); shading interp
plot([0,0],[freq(1)-0.25,freq(end)+0.25],'color',col.gray); axis tight
for h = find([clusStruct_catch1.clus_h]==1)
    dum = zeros(length(time_catch{1}),length(freq));
    [~,iat] = ismember(clusStruct_catch1(h).times,time_catch{1});
    [~,iaf] = ismember(clusStruct_catch1(h).freqs,freq);
    dum(sub2ind(size(dum), iat,iaf)) = 1;
    if clusStruct_catch1(h).peakZscore > 0
        cluscol = map(end,:);
    else
        cluscol = map(1,:);
    end
    contour(time_catch{1}, freq+[0.05,diff(freq)/2], dum',1, 'linewidth', 2, 'color', cluscol)
end
axis xy; caxis(cax); ylim(yax);
set(gca,'ytick','')
xlabel('Time (s)')
title('Catch - onset')
subplot(3,15,[13:14,28:29]); hold on;
pcolor(time_catch{2}, freq, squeeze(t_catch2)); shading interp
plot([0,0],[freq(1)-0.25,freq(end)+0.25],'color',col.gray); axis tight
for h = find([clusStruct_catch2.clus_h]==1)
    dum = zeros(length(time_catch{2}),length(freq));
    [~,iat] = ismember(clusStruct_catch2(h).times,time_catch{2});
    [~,iaf] = ismember(clusStruct_catch2(h).freqs,freq);
    dum(sub2ind(size(dum), iat,iaf)) = 1;
    if clusStruct_catch2(h).peakZscore > 0
        cluscol = map(end,:);
    else
        cluscol = map(1,:);
    end
    contour(time_catch{2}, freq+[0.05,diff(freq)/2], dum',1, 'linewidth', 2, 'color', cluscol)
end
axis xy; caxis(cax); ylim(yax);
set(gca,'ytick','')
xlabel('Time (s)')
title('Catch - response')

subplot(3,15,[15,30]); hold on;
imagesc(1, linspace(cax(1),cax(2),100), [1:100]')
set(gca,'xtick','','YAxisLocation','right')
ylabel('PPC_{correct} - PPC_{incorrect} (t-score)')
ylim(cax)

% plot RTs

xax = linspace(0,15,151);
dr = diff(xax(1:2));

subplot(3,15,31:32); hold on
plot([0,0],rax,'color',col.gray);
dumi = hist(RT_encoding_incorrect_cue,xax);
dumc = hist(RT_encoding_correct_cue,xax);
dumt = dumi + dumc;
plot(xax,dumi/sum(dumt)/dr, 'color',col.purple,'linewidth',2);
plot(xax,dumc/sum(dumt)/dr, 'color',col.green,'linewidth',2);
xlim([time_encoding{1}(1),time_encoding{1}(end)])
ylim(rax)
set(gca,'ytick',[0,0.5]);
xlabel('Time (s)')
ylabel('Response rate (Hz)')

subplot(3,15,33:34); hold on
plot([0,0],rax,'color',col.gray);
dumi = hist(RT_encoding_incorrect,xax);
dumc = hist(RT_encoding_correct,xax);
dumt = dumi + dumc;
plot(xax,dumi/sum(dumt)/dr, 'color',col.purple,'linewidth',2);
plot(xax,dumc/sum(dumt)/dr, 'color',col.green,'linewidth',2);
xlim([time_encoding{2}(1),time_encoding{2}(end)])
ylim(rax)
set(gca,'ytick',[]);
xlabel('Time (s)')

subplot(3,15,35:36); hold on
plot([0,0],rax,'color',col.gray);
plot(-xax,dumi/sum(dumt)/dr, 'color',col.purple,'linewidth',2);
plot(-xax,dumc/sum(dumt)/dr, 'color',col.green,'linewidth',2);
xlim([time_encoding{3}(1),time_encoding{3}(end)])
ylim(rax);
set(gca,'ytick',[],'yticklabel',{});
xlabel('Time (s)')

subplot(3,15,37:38); hold on
plot([0,0],rax,'color',col.gray);
dumi = hist(RT_retrieval_incorrect,xax);
dumc = hist(RT_retrieval_correct,xax);
dumt = dumi + dumc;
plot(xax,dumi/sum(dumt)/dr, 'color',col.purple ,'linewidth',2);
plot(xax,dumc/sum(dumt)/dr, 'color',col.green ,'linewidth',2);
xlim([time_retrieval{1}(1),time_retrieval{1}(end)])
ylim(rax)
set(gca,'ytick',[],'yticklabel',{});
xlabel('Time (s)')

subplot(3,15,39:40); hold on
plot([0,0],rax,'color',col.gray);
plot(-xax,dumi/sum(dumt)/dr, 'color',col.purple,'linewidth',2);
plot(-xax,dumc/sum(dumt)/dr, 'color',col.green,'linewidth',2);
xlim([time_retrieval{2}(1),time_retrieval{2}(end)])
ylim(rax)
set(gca,'ytick',[],'yticklabel',{});
xlabel('Time (s)')

subplot(3,15,41:42); hold on
plot([0,0],rax,'color',col.gray);
dumi = hist(RT_catch_incorrect,xax);
dumc = hist(RT_catch_correct,xax);
plot(xax,dumi/sum(dumt)/dr, 'color',col.purple,'linewidth',2);
plot(xax,dumc/sum(dumt)/dr, 'color',col.green,'linewidth',2);
xlim([time_catch{1}(1),time_catch{1}(end)])
ylim(rax)
set(gca,'ytick',[],'yticklabel',{});
xlabel('Time (s)')

subplot(3,15,43:44); hold on
plot([0,0],rax,'color',col.gray);
plot(-xax,dumi/sum(dumt)/dr, 'color',col.purple,'linewidth',2);
plot(-xax,dumc/sum(dumt)/dr, 'color',col.green,'linewidth',2);
xlim([time_catch{2}(1),time_catch{2}(end)])
ylim(rax)
set(gca,'ytick',[],'yticklabel',{});
xlabel('Time (s)')
