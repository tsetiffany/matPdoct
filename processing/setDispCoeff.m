function dispCoeff = setDispCoeff(inputData, depths, maxDispOrders, coeffRange)

arrCountDispCoeff = zeros(1, maxDispOrders-1);

frameData = inputData;

for i = 1 : length(arrCountDispCoeff)
%     arrDispCoeffRng = [coeffRange, -1 * coeffRange];
    arrDispCoeffRng = linspace(-1 * coeffRange,coeffRange,50);
    arrCost = zeros(size(arrDispCoeffRng));

    for j = 1 : length(arrDispCoeffRng)
        arrCountDispCoeff(i) = arrDispCoeffRng(j);
        arrCost(j) = calCostFunc(frameData, depths, maxDispOrders, arrCountDispCoeff);
    end
    
    for k = 1 : 50
        [dispCoeff_min1, dispCoeff_min2] = srcMinCostCoeff(arrDispCoeffRng, arrCost);
        dispCoeff_new = (dispCoeff_min1 + dispCoeff_min2)/2;
        arrDispCoeffRng = [arrDispCoeffRng, dispCoeff_new];
        arrCountDispCoeff(i) = dispCoeff_new;
        cost_new = calCostFunc(frameData, depths, maxDispOrders, arrCountDispCoeff);
        arrCost = [arrCost, cost_new];
    end
    [~, argmin] = min(arrCost);
    arrCountDispCoeff(i) = arrDispCoeffRng(argmin);
end

dispCoeff = arrCountDispCoeff;
end
