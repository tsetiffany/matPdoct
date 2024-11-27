function RawData_PhaseComp = compPhase(fftData_2D, phaseSlope)

pixIdx = linspace(1,size(fftData_2D,1),size(fftData_2D,1));

for i = 1:size(fftData_2D, 2)
    ref_FFTData_Comp(:,i) = fftData_2D(:,i)...
        .*exp(-1j*(pixIdx'*phaseSlope(i)));
end

RawData_PhaseComp = ifft(ref_FFTData_Comp);

end