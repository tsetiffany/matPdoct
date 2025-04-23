function [frstMin, scndMin] = srcMinCostCoeff(arrCoeff, arrCost)

% arrCoeff : array of dispersion coefficient
% arrCost : array of cost value
% frstMin : first minimum cost value
% scndMin : second minimum cost value

[val, idx] = sort(arrCost);
frstMin = arrCoeff(idx(1));
scndMin = arrCoeff(idx(2));

end