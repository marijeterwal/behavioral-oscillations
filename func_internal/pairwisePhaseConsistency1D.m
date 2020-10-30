function ppc = pairwisePhaseConsistency1D(data, dim)
%{ 
PairwisePhaseConsistency1D

Computes Pairwise Phase Consistency (PPC) on angle data along one dimension dim.

For more information on Pairwise Phase Consistency see:
Vinck, M., van Wingerden, M., Womelsdorf, T., Fries, P., and Pennartz,
C.M.A. (2010). The pairwise phase consistency: a bias-free measure of 
rhythmic neuronal synchronization. Neuroimage 51, 112–122.

Marije ter Wal - 2020
m.j.terwal@bham.ac.uk

%}

if isempty(dim)
    dim = 1;
end

dims = size(data);
dims(dim) = 1;

dataPrep = permute(data,[dim,setdiff(1:length(dims),dim)]);

if length(dim) > 1
    error('More than one dimension specified; this option is not supported.')
elseif length(dim) == 1
    N = zeros(size(dataPrep(1,:)));
    dumsum = zeros(size(dataPrep(1,:)));
    for p1 = 1:size(data,dim)-1
        for p2 = p1+1:size(data,dim)    
        dum = cos(dataPrep(p1,:) - dataPrep(p2,:));
            
            N(~isnan(dum)) = N(~isnan(dum)) + 1;
            dum(isnan(dum)) = 0;
            
            dumsum = dumsum + dum;
        end
    end
    ppc = reshape(dumsum./N, dims);
end

end