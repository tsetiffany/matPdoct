function phaseSlope = estPhaseShift(ref_fftData_1D, fftData_2D, offSetIdx)
% ref_fftData_1D = ref_Ascan_1;
% fftData_2D = ref_FFTData_rescaled(:,1:4:end);
% offSetIdx = calSigOffSetIdx;

cplxConjCalSig  = fftData_2D.*repmat(conj(ref_fftData_1D), [1 size(fftData_2D,2)]);
[~, idx]        = sort(abs(ref_fftData_1D(offSetIdx:round(end/2))), 'descend');
CalSigIdx       = offSetIdx+idx(1)-1;
diffPhases      = angle(cplxConjCalSig(CalSigIdx,:));

pixIdx = linspace(1,size(ref_fftData_1D,1),size(ref_fftData_1D,1));

phaseSlope = diffPhases /CalSigIdx;

end