function gPhaseOffset = calGlobPhaseOff(JM11, JM12, JM21, JM22)
%calGlobPhaseOff(JM11, JM12, JM21, JM22) - calculates global phase offset
%
% The function takes the four entries of the JM in 3D matrix form (JM##,
% size=[rFrames,width,depth]). Returns the global phase offset
% (gPhaseOffset, size=[rFrames,width,depth])
%
% Using rFrame=1 as the reference JM matrix, the amplitude and phase for
% each entry is computed separately, then combined into polar complex forms
% for each gpJM##.  The final gPhaseOffset is then found by summing the
% complex gpJM## values and taking the angle.
%

% get the delta phi(o,j) amplitudes w.r.t. the first rFrame
gpJM11_amp  = 1/((1/abs(JM11)) + (1/abs(JM11(1,:,:))));
gpJM12_amp  = 1/((1/abs(JM12)) + (1/abs(JM12(1,:,:))));
gpJM21_amp  = 1/((1/abs(JM21)) + (1/abs(JM21(1,:,:))));
gpJM22_amp  = 1/((1/abs(JM22)) + (1/abs(JM22(1,:,:))));

% get the delta phi(o,j)phases w.r.t. the first rFrame
gpJM11_phs  = angle(JM11) - angle(JM11(1,:,:));
gpJM12_phs  = angle(JM12) - angle(JM12(1,:,:));
gpJM21_phs  = angle(JM21) - angle(JM21(1,:,:));
gpJM22_phs  = angle(JM22) - angle(JM22(1,:,:));


% put the complex values in polar (|amp|*exp(iPhs)) format
gpJM11_cplx = gpJM11_amp.*exp(i.*gpJM11_phs);
gpJM12_cplx = gpJM12_amp.*exp(i.*gpJM12_phs);
gpJM21_cplx = gpJM21_amp.*exp(i.*gpJM21_phs);
gpJM22_cplx = gpJM22_amp.*exp(i.*gpJM22_phs);

% get the angle of their sum
gPhaseOffset = angle(gpJM11_cplx+gpJM12_cplx+gpJM21_cplx+gpJM22_cplx);

end