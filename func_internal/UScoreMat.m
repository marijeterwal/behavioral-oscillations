function [Z,sgnf] = UScoreMat(dat, ref, alpha)
%{ 
----- UscoreMat -----

DESCRIPTION:
Computes Mann-Whitney U-test between data matrices dat and ref and outputs
a Z-score approximation and significance at level alpha. Every data point
in dat is compared to the reference data in the rows of ref.

INPUTS:
- dat: The function expects a data matrix dat of size n1 x n2 x n3 x...
- ref: The first dimension of ref is expected to contain repetitions against which 
the data in dat have to be compared (i.e. nrep x n1 x n2 x n3 x ...)
- alpha: significance threshold

OUTPUTS:
- Z: Z-scored data, same size as dat
- sgnf: significance at level alpha (0 or 1)

Marije ter Wal - 2020
m.j.terwal@bham.ac.uk

%}

dat = reshape(dat, [1,size(dat,1), size(dat,2)]);

%% normal Mann Whitney U test

n2 = size(ref,1);
n = n2+1;

U = squeeze(sum(bsxfun(@lt,sort(ref,1),dat),1) + 1);

mU = n2/2;
sigU = sqrt(n2*(n2+2)/12);

Z = (U - mU) / sigU; 

%% normal Z-test

sgnf = sign(Z) .* (abs(Z) >= norminv(1-alpha, 0, 1));

end