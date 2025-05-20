function gPhaseOffset = calGlobPhaseOff_2JMs(JM11, JM12)
%calGlobPhaseOff(JM11, JM12) - calculates global phase offset
%
% The function takes the four entries of the JM in 3D matrix form (JM##,
% size=[depth, width, rFrames]). Returns the global phase offset
% (gPhaseOffset, size=[depth, width, rFrames])
%
% Using rFrame=1 as the reference JM matrix, the amplitude and phase for
% each entry is computed separately, then combined into polar complex forms
% for each gpJM##.  The final gPhaseOffset is then found by summing the
% complex gpJM## values and taking the angle.
%

% % iterate over each rFrame in JM##
% gpJM11_amp = zeros(size(JM11));
% gpJM12_amp = zeros(size(JM12));
% gpJM11_phs = zeros(size(JM11));
% gpJM12_phs = zeros(size(JM12));
% 
% for idx = 1:size(JM11, 3)
%     % get the delta phi(o,j) amplitudes w.r.t. the first rFrame
%     gpJM11_amp(:,:,idx)  = 1./((1./abs(JM11(:,:,idx))) + (1./abs(JM11(:,:,1))));
%     gpJM12_amp(:,:,idx)  = 1./((1./abs(JM12(:,:,idx))) + (1./abs(JM12(:,:,1))));
% 
%     % get the delta phi(o,j)phases w.r.t. the first rFrame
%     gpJM11_phs(:,:,idx)  = angle(JM11(:,:,idx)) - angle(JM11(:,:,1));
%     gpJM12_phs(:,:,idx)  = angle(JM12(:,:,idx)) - angle(JM12(:,:,1));
% 
% end

% get the delta phi(o,j) amplitudes w.r.t. the first rFrame
gpJM11_amp  = 1./((1./abs(JM11)) + (1./abs(JM11(:,:,1))));
gpJM12_amp  = 1./((1./abs(JM12)) + (1./abs(JM12(:,:,1))));


% get the delta phi(o,j)phases w.r.t. the first rFrame
gpJM11_phs  = angle(JM11) - angle(JM11(:,:,1));
gpJM12_phs  = angle(JM12) - angle(JM12(:,:,1));


% put the complex values in polar (|amp|*exp(iPhs)) format
gpJM11_cplx = gpJM11_amp.*exp(1i.*gpJM11_phs);
gpJM12_cplx = gpJM12_amp.*exp(1i.*gpJM12_phs);


% get the angle of their sum
gPhaseOffset = angle(gpJM11_cplx+gpJM12_cplx);

end