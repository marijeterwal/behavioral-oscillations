function [out] = extractClusters(data, thres)
%{ 
----- extractClusters -----

This function uses the Connected-Component Labeling approach (Dillencourt 
et al., 1992) to identify unique connected areas in data above threshold thres.

INPUTS:
- data is a 2D data array
- thres is the binarization threshold

OUTPUT:
- 2D array of the size of data specifying the cluster labels (0 means data
was below threshold)

DESCRIPTION:
The data in dat are thresholded at thres.
Then, the 2-pass algorithm walks through the pixel space and compares all 
pixels to the values of its four top and left neighbors (‘8-connectivity’). 
In the first pass, each pixel with value 1 or -1 receives either a unique 
label, if the pixel has no neighbors with the same value, or with the same 
label as its left neighbor. In the second pass, connected regions with 
different labels are merged, leaving unique regions with unique labels.

For details see: 
Dillencourt, M. B., Samet, H., & Tamminen, M. (1992). A general approach to 
connected-component labeling for arbitrary image representations. Journal 
of the ACM, 39(2), 253–280. http://doi.org/10.1145/128749.128750

Marije ter Wal - 2020
m.j.terwal@bham.ac.uk

%}

if ~isempty(thres)
    % binarize image
    bindata = zeros(size(data));
    bindata(data >= thres) = 1;%2;
    % bindata(data <= -thres) = 1;
else
    bindata = data;
end

% set boundaries
bindata(1,:) = 0;
bindata(end,:) = 0;
bindata(:,1) = 0;
bindata(:,end) = 0;
 
nextlabel = 3;
linked = {}; linked{1} = [1]; linked{2} = [2]; linked{3} = [3];

% first pass
for rw = 2:size(data,1)-1
    for cl = 2:size(data,2)-1
        
        if bindata(rw,cl) ~= 0
            % find neighbors
            neighbors = [bindata(rw-1,cl-1:cl+1), bindata(rw,cl-1)];
            neighborslb = neighbors(find(neighbors));
            
            % are neighbors belonging to the same group?
            if isempty(find(neighborslb ~= bindata(rw,cl))) % no - create new group
                bindata(rw,cl) = nextlabel;
                linked{nextlabel} = [nextlabel];
                nextlabel = nextlabel + 1;
            else % yes - link with neighbors
                % find smallest label
                bindata(rw,cl) = nanmin(neighborslb);
                for ul = unique(neighborslb)
                    linked{ul} = unique(horzcat(linked{ul}, unique(neighborslb), linked{min(unique(neighborslb))}));
                end
            end  
        end
    end
end

% second pass - merge labels that were found to be identical
for rw = 2:size(data,1)-1
    for cl = 2:size(data,2)-1
        if bindata(rw,cl) ~= 0
            bindata(rw,cl) = min(linked{bindata(rw,cl)});
        end
    end
end

out = bindata;
end