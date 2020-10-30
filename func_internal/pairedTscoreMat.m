function [t,nu] = pairedTscoreMat(dat1,dat2,dim)
%{ 
----- pairedTscoreMat -----

DESCRIPTION:
Computes the paired t-score of data dat1 and dat2 along dimension dim

INPUTS:
- dat1: data array
- dat2: data array (optional). If not specified, function computes
one-sample t-scoring
- dim: dimensional along which data need to be t-scored.

OUTPUTS:
- t: t-scores (array size as in dat1)
- nu: degrees of freedom (array size as in dat1)


Marije ter Wal - 2020
m.j.terwal@bham.ac.uk

%}


if ~isempty(dat2)
    assert(sum(size(dat1) == size(dat2)) == length(size(dat1)))
    dat = dat1-dat2;
else
    dat = dat1;
end
n = size(dat,dim);

t = nanmean(dat,dim) ./ (nanstd(dat,0,dim)/sqrt(n));
nu = (n-1) * ones(size(t));

end