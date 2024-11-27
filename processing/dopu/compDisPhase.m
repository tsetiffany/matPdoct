function frameData_dispComp = compDisPhase(frameData, maxDispOrders, arrCountDispCoeff)

ScanPts         = size(frameData, 1);
LinePerFrame    = size(frameData, 2);
kLinear         = linspace(-1,1,ScanPts);
kaxis           = repmat(kLinear',1,LinePerFrame);
amp             = abs(frameData);
phase           = angle(frameData);

for i = 1:maxDispOrders-1
    phase = phase + arrCountDispCoeff(i)*(kaxis.^(i+1));

end

frameData_dispComp = amp.*exp(1j.*phase);

end
