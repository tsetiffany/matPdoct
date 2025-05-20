function cost = calCostFun(frameData, depths, maxDispOrders, arrCountDispCoeff)

frameData_dispComp = compDisPhase(frameData,maxDispOrders,arrCountDispCoeff);

OCT      = ((abs(fft(frameData_dispComp))).^2);
roiOCT   = OCT(depths(1):depths(2),:);
% imagesc(imadjust(mat2gray(roiOCT))); colormap(gray)
% pause(0.01)
normOCT  = roiOCT./sum(roiOCT(:));
entropy  = -1*((normOCT).*log(normOCT));
cost     = sum(entropy(:));
end




