function clusStruct = detectClusters(cfg, data1, threshold, data2, dataRef)
%{ 
----- detectClusters -----

This function identifies clusters in 2D data and splits them into 
individual peaks.

INPUTS:
- cfg: struct with settings (see below)
- data is a 2D data array
- thres is the binarization threshold
- dataZ (optional)
- dataRef (optional)

OUTPUT:
- clusStruct: a struct with, for every cluster, a series of information,
such as the times and frequencies of all pixels in the cluster, the center
of mass of the cluster, the peak and sum of the cluster values, etc. 

DESCRIPTION:
The first step in this function is the extraction of clusters in the data
that are above or below threshold. The function uses one of two methods for
extracting the clusters: 
- 'extract': based on an 8-connectivity connected component labeling
approach
- 'watershed': based on a watershed algorithm

If using watershed, the data should be smoothed first, for the extract
method this is optional.

The function then goes through the following optional steps:
- split the clusters into individual peaks
- take and store the cluster data from a separate dataset data2 instead of
from data1
- store info about data1/data2 as well as about the same cluster locations 
for a separate reference dataset dataRef 

Splitting of clusters is done through the following steps:
- the threshold is iteratively increased towards the highest value in the 
clusters;
- the cluster extraction is redone for every threshold value;
- The number of subclusters is assessed. Resulting subclusters are 
required to be at least the fraction cfg.minClusterSize of the original cluster
to be considered cluster (this overcomes noise in the data).
- If no subclusters of sufficient size are found, the threshold is increased 
further. If all identified subclusters were smaller than cfg.minClusterSize
cluster splitting is aborted. 
- If subclusters of sufficient size are present, these are stored. 
- All pixels that were part of the original cluster, but are not a member 
of any of the new subclusters, the weighted Euclidian distance to all new 
subclusters is computed. The pixels are assigned to the closest subcluster. 

CONFIG:
- cfg.preEvent: start of the time series
- cfg.postEvent: last time point in time series
- cfg.freqs: frequencies (corresponding to the rows in data)
- cfg.method: cluster extraction method: 'extract' or 'watershed'
- cfg.smooth: double or array defining the smoothing width of the gaussian
filter (in indices)
- cfg.splitClusters: boolean to determine whether clusters are split into
peaks
- cfg.peakHeight: (required for splitClusters) vector of relative peak 
heights used as threshold in cluster splitting
- cfg.minNewClusSize: (required for splitClusters) minimum cluster size of 
split clusters


Marije ter Wal - 2020
m.j.terwal@bham.ac.uk

%}


% prepare output structure
clusStruct = struct('clusID',[], 'times', [],...
    'freqs',[], 'CoMtime', [], 'CoMfreq',[], 'pAdj',[], 'peakZscore', [], 'peakpAdj',[]);

% check config
tspan       = cfg.preEvent:cfg.dt:cfg.postEvent;
if isfield(cfg,'freqs')
    freqs = cfg.freqs;
elseif isfield(cfg,'scal') && isfield(cfg,'cw') % for backwards compatibility
    wav     = ['cmor', num2str(cfg.cw), '-1'];
    scal    = cfg.scal;
    freqs   = scal2frq(scal, wav, 1/cfg.fs);
else
    error('cfg is incomplete') % elaborate
end

if ~isfield(cfg,'splitClusters')
    cfg.splitClusters = false;
else
    if ~isfield(cfg,'peakHeight')
        cfg.peakHeight = 0.05:0.05:0.5; % pixels
    end
    if ~isfield(cfg,'minNewClusSize')
        cfg.minNewClusSize = 0.1;
    end
end
if ~isfield(cfg,'method') || isempty(cfg.method)
    cfg.method = 'watershed';
end
if ~isfield(cfg,'smooth')
    cfg.smooth = [];
elseif cfg.smooth == true
    cfg.smooth = [1,1];
end

if strcmp(cfg.method, 'extract')
    % smooth
    if cfg.smooth
        sdata = imgaussfilt(data1,cfg.smooth);
    else
        sdata = data1;
    end
    
    % cluster extraction
    clusDatPlus = extractClusters(sdata, threshold);
    clusDatMin = extractClusters(-1*sdata, threshold);
    clusDat = -1*clusDatMin + clusDatPlus;
    
elseif strcmp(cfg.method,'watershed')
    
    % ----- positive clusters
    % threshold image 
    tdata = -1*data1;
    tdata(data1<threshold) = 0;

    % smooth image
    sdata = imgaussfilt(tdata,cfg.smooth);

    % find labels
    clusDatPlus = watershed(100*sdata,8);
    clusDatPlus(~sdata) = 0;
    
    % ----- negative clusters
    tdata = data1;
    tdata(-1*data1>threshold) = 0;
    
    % smooth image
    sdata = imgaussfilt(tdata,cfg.smooth);
    
    % find labels
    clusDatMin = watershed(sdata,8);
    clusDatMin(~sdata) = 0;
    
    % ----- all clusters
    clusDat = -1*clusDatMin + clusDatPlus;
end

%% split clusters with multiple peaks

if isfield(cfg,'splitClusters') && cfg.splitClusters == true
    
    clusLabels = setdiff(unique(clusDat),0);
    
    maxlabelPlus = max(clusDat(clusDat(:)>=0));
    maxlabelMin = min(clusDat(clusDat(:)<0));
    
    % look for peaks within each cluster
    bi = 1;
    while bi <= length(clusLabels)
        %     for bi = 1:length(clusLabels)
        
        if clusLabels(bi)>0
            fact = 1;
        else
            fact = -1;
        end
        
        cont = true;
        cnt = 1;
        while cont && cnt<=length(cfg.peakHeight)
            mm = max(abs(data1(clusDat==clusLabels(bi))));
            dataRed = fact*sdata;
            dataRed = dataRed - mm*cfg.peakHeight(cnt);
            dataRed(clusDat~=clusLabels(bi)) = 0;
            
            if sum(dataRed(:)) == 0
                cont = false; % cluster of size 1
            else
                clusField = extractClusters(dataRed,threshold);
                newClusLabels = setdiff(unique(clusField),0);
                
                % to determine whether this setting is correct
                newClusSize = arrayfun(@(x) sum(clusField(:)==x), newClusLabels);
                oldClusSize = sum(clusDat(:)==clusLabels(bi));
                
                if length(newClusLabels) < 2
                    % nothing changed, move on to next setting
                    %                     cont = false;
                    cnt = cnt+1;
                    continue
                elseif sum(newClusSize/oldClusSize > cfg.minNewClusSize) == 0
                    % all new clusters are too small, move on to next
                    % cluster
                    cont = false;
                elseif sum(newClusSize/oldClusSize > cfg.minNewClusSize) == 1
                    % this cluster setting is not appropriate, try one step
                    % higher
                    cnt = cnt+1;
                    continue
                else
                    % stick with this peakHeight setting and update clusters
                    
                    % take all new clusters larger than minNewClusSize
                    newClusLabels = newClusLabels(newClusSize/oldClusSize > cfg.minNewClusSize);
                    clusField(~ismember(clusField,newClusLabels)) = 0;
                    
                     dataFull = clusDat==clusLabels(bi);
                     dataMissing = dataFull - clusField>0;
                     
                     [missIDy, missIDx] = ind2sub(size(dataMissing),find(dataMissing));
                     if ~isempty(missIDy)
                         
                     i2 = 1;
                     dist = [];
                     for bii = newClusLabels'
                         [fieldIDy, fieldIDx] = ind2sub(size(clusField),find(clusField==bii));
                         fp = []; tp = [];
                         for i = 1:length(missIDy)
                             fp(i,:) = (freqs(missIDy(i))- freqs(fieldIDy)) / (freqs(end)-freqs(1));
                             tp(i,:) = (tspan(missIDx(i))- tspan(fieldIDx)) / (tspan(end)-tspan(1));
                         end
                         % euclidean distance to the cluster
                         dist(i2,:) = min(sqrt(cfg.freqweight*fp.^2 + cfg.timeweight*tp.^2),[],2); %  (length(tspan)/length(freqs)) clusters x pixels
                         i2 = i2 + 1;
                     end
                     [mindist,id] = min(dist,[],1);
                     clusField(sub2ind(size(clusField),missIDy,missIDx)) = newClusLabels(id);
                     end
                 
                    % update labels
                    if clusLabels(bi) > 0
                        clusDat(clusDat==clusLabels(bi)) = clusDat(clusDat==clusLabels(bi)) + maxlabelPlus + clusField(clusDat==clusLabels(bi));
                        maxlabelPlus = max(clusDat(clusDat(:)>=0));
                    else
                        clusDat(clusDat==clusLabels(bi)) = clusDat(clusDat==clusLabels(bi)) + maxlabelMin - clusField(clusDat==clusLabels(bi));
                        maxlabelMin = min(clusDat(clusDat(:)<0));
                    end
                    cont = false;
                end
            end
        end
        bi = bi+1;
        clusLabels = setdiff(unique(clusDat),0);
    end
end

%% characterize each detected blob...

clusLabels = setdiff(unique(clusDat),0);

if length(clusLabels) == 0
    clusStruct = [];
    return
end

for bi = 1:length(clusLabels)
    clusid = clusLabels(bi);
    cluslocs = find(clusDat == clusid);
    [I,J] = ind2sub(size(clusDat),cluslocs);
    com_scale = sum(I)/length(I);
    com_time = sum(J)/length(J);
  
    if ~isempty(data2)
        [~, peakid] = max(abs(data2(cluslocs)));
        peakZ = data2(cluslocs(peakid));
        pAdj = FDRfromZscore(data2, cfg.qval); % compute blob FDR
        Zclus = data2(cluslocs);
    else
        [~, peakid] = max(abs(data1(cluslocs)));
        peakZ = data1(cluslocs(peakid));
        sumZ = sum(sum(data1(cluslocs)));
        if isfield(cfg, 'qval')
            pAdj = FDRfromZscore(data1, cfg.qval);
        else
            pAdj = NaN(size(data1));
        end
        Zclus = data1(cluslocs);
    end
    
    % save everything to struct
    clusStruct(bi).blobID = clusid;
    clusStruct(bi).times = tspan(J);
    clusStruct(bi).freqs = freqs(I);
    clusStruct(bi).Zscores = Zclus;
    clusStruct(bi).CoMtime = tspan(round(com_time));
    clusStruct(bi).CoMfreq = freqs(floor(com_scale));
    clusStruct(bi).pAdj = pAdj(cluslocs);
    clusStruct(bi).peakZscore = peakZ;
    clusStruct(bi).sumZscore = sumZ;
    clusStruct(bi).peakpAdj = pAdj(cluslocs(peakid));
    
    if ~isempty(dataRef)
        clusStruct(bi).peakZscoreRef = dataRef(cluslocs(peakid));
        clusStruct(bi).maxZscoreRef = nanmax(dataRef(cluslocs));
        clusStruct(bi).avgZscoreRef = nanmean(dataRef(cluslocs));
    end
    
end

end

function [adj_p,h] = FDRfromZscore(data, qval)
% data is a set (n-D vector) of data points

% find p-value
p = normcdf(-abs(data));

[tmp1, ~, ~, tmp2] = fdr_bh(p(:), qval,'pdep');

adj_p = reshape(tmp2,size(data));
h = reshape(tmp1,size(data));

end