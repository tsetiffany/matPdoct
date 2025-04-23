function frameData_dispComp = compDisPhase(frameData, maxDispOrders, arrCountDispCoeff)

ScanPts         = size(frameData, 1);
LinePerFrame    = size(frameData, 2);
kLinear         = linspace(-1,1,ScanPts);
kaxis           = repmat(kLinear',1,LinePerFrame);

frameData_dispComp = frameData;
for i = 1:maxDispOrders-1
    frameData_dispComp = frameData_dispComp.*exp(1j.*(arrCountDispCoeff(i)*(kaxis.^(i+1))));
end

end

