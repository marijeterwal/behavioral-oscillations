function [h,p] = nonParamPVal(data,ref,alpha)
%{ 
----- nonParamPVal -----

DESCRIPTION:
Computes the fraction of data points in reference data set ref that are
larger (or smaller) than the original data. 

INPUTS:
- data: data array
- ref: reference data array
- alpha: significance level

OUTPUTS:
- h: significance at level alpha
- p: p-value


Marije ter Wal - 2020
m.j.terwal@bham.ac.uk

%}


% sort ref
[vals, ~] = sort(ref); 

p = zeros(length(data),1);
h = zeros(length(data),1);

for c = 1:length(data)
    mu = median(vals);
    if data(c) < mu
        dum = find(data(c) <= vals, 1, 'first')-1;
    else
        dum = length(vals) - find(data(c) >= vals, 1, 'last') -1;
    end
    p(c) = dum/length(vals);
    
    if p(c)<= alpha
        h(c) = 1;
    else
        h(c) = 0;
    end
end

end