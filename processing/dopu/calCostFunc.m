function cost = calCostFun(frameData, maxDispOrders, arrCountDispCoeff)

frameData_dispComp = compDisPhase(frameData,maxDispOrders,arrCountDispCoeff);

OCT      = abs(fft(frameData_dispComp)).^2;
roiOCT   = OCT(150:350,:);
normOCT  = roiOCT./sum(roiOCT(:));
entropy  = -1*((normOCT).*log(normOCT));
cost     = sum(entropy(:));
end




