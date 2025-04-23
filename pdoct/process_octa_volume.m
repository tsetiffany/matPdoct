function [avgOCT_mcorr, avgOCT_tcorr, OCTA_mcorr, OCTA_tcorr] = process_octa_volume(cplxData_A,cplxData_B,numMscans,process_path)

% Location of motion correction code
addpath('C:\Users\tiffa\Documents\1. Projects\PDOCT Processing\main_vctrl\Code_workingversion\matPdoct\MotionCorrection');

[numPoints,numAlines,numBscans] = size(cplxData_A);

% Bulk-phase correction between P and S channels
for i = 1:numBscans
    OCT_P  = cplxData_A(:,:,i);
    OCT_S  = cplxData_B(:,:,i);

    rPhaseOff = repmat(angle(sum(OCT_S.*conj(OCT_P),1)), [size(OCT_S,1) 1]);
    OCT_Cplx = (OCT_P) + (OCT_S.*exp(-1j.*rPhaseOff));
    CompCplx(:,:,i) = OCT_Cplx;
end

save(fullfile(process_path,'CompCplx.mat'), 'CompCplx', '-v7.3');

clearvars cplxData_A cplxData_B

% Remove 1st BM scan in each set
frames_to_keep  = setdiff(1:size(CompCplx, 3), 1:4:size(CompCplx, 3)); % for 4BM
CompCplx_BM     = CompCplx(:, :, frames_to_keep);
CompCplx_BM_ROI = CompCplx_BM(:, 26:end-25,:); % cut flyback

% Local sub-pixel motion correction 
disp('start sub-pixel motion correction')
numMscans = 3; % set to 3 after removing 1st BM scan
cplxOCT_mcorr  = CompCplx_BM_ROI;
numFrames = size(CompCplx_BM_ROI, 3);
for I = 1:numMscans:numFrames-1
    batchFrames = cplxOCT_mcorr(:, :, I:I+numMscans-1);
    correctedBatch = MotionCorrection(batchFrames);
    cplxOCT_mcorr(:, :, I:I+numMscans-1) = correctedBatch;
end

disp('finished sub-pixel motion correction')
clearvars -except cplxOCT_mcorr 

% Bulk phase offset between BM scans
numMscans = 3;
iF=1; % frame counter
for I = 1:numMscans:size(cplxOCT_mcorr,3)
    K = ((I-1)/numMscans)+1;
    BMscan     = cplxOCT_mcorr(:,:,I:I+(numMscans-1));
    BMscan_sub = cplxOCT_mcorr(:,:,I:I+(numMscans-1)-1);
    for J = 1:numMscans-1
        Xconj   = BMscan(:,:,1+J).*conj(BMscan(:,:,1));
        BulkOff = repmat(angle(sum(Xconj,1)), [size(Xconj,1) 1]);
        BMscan(:,:,1+J)   = BMscan(:,:,1+J).*exp(-1j*BulkOff);
        BMscan_sub(:,:,J) = BMscan(:,:,1) - BMscan(:,:,1+J);
    end
    % Average OCT %
    avgOCT(:,:,K) = abs(mean(BMscan,3));
    % Variance %
%     VAR(:,:,K) = abs(var(BMscan,0,3));
    % Subtraction %
    SUB(:,:,K) = mean(abs(BMscan_sub),3);
    % Complex Differential Variance %
%     CDV(:,:,K) = abs(var(BMscan_sub,0,3));
    if mod(iF,10) == 0
        fprintf('OCTA main process : %d\n', iF);
    end
    iF = iF + 1;
end


% % Local sub-pixel motion correction 
% disp('start sub-pixel motion correction')
% numMscans = 3; % set to 3 after removing 1st BM scan
% cplxOCT_mcorr  = CompCplx;
% numFrames = size(CompCplx, 3);
% for I = 1:numMscans:numFrames-1
%     batchFrames = cplxOCT_mcorr(:, :, I:I+numMscans-1);
%     correctedBatch = MotionCorrection(batchFrames);
%     cplxOCT_mcorr(:, :, I:I+numMscans-1) = correctedBatch;
% end
% 
% disp('finished sub-pixel motion correction')
% clearvars -except cplxOCT_mcorr 

% % Bulk phase offset between BM scans
% numMscans = 3;
% iF=1; % frame counter
% for I = 1:numMscans:750%size(cplxOCT_mcorr,3)
%     K = ((I-1)/numMscans)+1;
% 
%     Xconj_2   = cplxOCT_mcorr(:,:,I+1).*conj(cplxOCT_mcorr(:,:,I));
%     Xconj_3   = cplxOCT_mcorr(:,:,I+2).*conj(cplxOCT_mcorr(:,:,I));
% 
%     BulkOff_2 = repmat(angle(sum(Xconj_2)), [size(Xconj_2,1) 1]);
%     BulkOff_3 = repmat(angle(sum(Xconj_3)), [size(Xconj_3,1) 1]);
% 
%     Bscan_1  = cplxOCT_mcorr(:,:,I);
%     Bscan_2  = cplxOCT_mcorr(:,:,I+1) .* exp(-1j*BulkOff_2);
%     Bscan_3  = cplxOCT_mcorr(:,:,I+2) .* exp(-1j*BulkOff_3);
% 
%     % Average OCT %
%     avgOCT(:,:,K) = (abs(Bscan_1) + abs(Bscan_2) + abs(Bscan_3))/3;
% 
% 
%     % Complex Subtraction
%     % 		Var(:,:,K) = abs(var(cat(3,(Bscan_1),(Bscan_2),(Bscan_3)),0,3));
%     Sub(:,:,K) = abs(Bscan_1 - Bscan_2) + abs(Bscan_2 - Bscan_3);
%     
%     sub1 = Bscan_1 - Bscan_2;
%     sub2 = Bscan_2 - Bscan_3;
% 
%     if mod(iF,10) == 0
%         fprintf('OCTA main process : %d\n', iF);
%     end
%     iF = iF + 1;
% end

clearvars -except avgOCT SUB

% Global Motion Correction
numFrames   = size(avgOCT,3);
yShift_axial = zeros([numFrames 1]);
limit = 150;
usfac = 1; 
avgOCT_mcorr = avgOCT;
OCTA_mcorr = SUB;

for I = 1:numFrames
    % Every 'for' loop, reference frame will be the middle frame
    [output, ~] = dftregistration_vol(fft2(20.*log10(avgOCT(:, :, round(numFrames./2)))), ...
        fft2(20.*log10(avgOCT(:, :, I))), usfac);

        if isempty(output)
            output = zeros(1,4);
        end

    % Assign and save the shifting value for axial (yShift)
    yShift_axial(I) = round(output(3));
    
    % Thresholding  value was found via plotting the shifting values
    if abs(output(3)) >= limit
        output(3) = 0;
    end

    avgOCT_mcorr(:,:,I)  = circshift(avgOCT(:,:,I), [yShift_axial(I), 0]);
    OCTA_mcorr(:,:,I)  = circshift(SUB(:,:,I), [yShift_axial(I), 0]);
end

% Tilt correction
numLines = size(avgOCT_mcorr,2);
usfac = 1;

avgOCT_tcorr = zeros(size(avgOCT_mcorr));
OCTA_tcorr = zeros(size(OCTA_mcorr));

for I = 1:numLines

    %Every 'for' loop, reference frame will be the middle frame
    [output,~] = dftregistration_vol(fft2(squeeze(avgOCT_mcorr(:, round(numLines./2), :))), fft2(squeeze(avgOCT_mcorr(:,I,:))), usfac);
    
    %Assign and save the shifting value along axial direction
    axialShift_tilt(I) = round(output(3));

end

%Curve Fitting
x = [1:numLines]';
cx = polyfit(x,axialShift_tilt(:,1:end),2);
axialShift_tilt_fit = polyval(cx,x);

for j = 1:numLines
    avgOCT_tcorr(:,j,:) = circshift(avgOCT_mcorr(:,j,:), [round(axialShift_tilt_fit(j)), 0]);
    OCTA_tcorr(:,j,:) = circshift(OCTA_mcorr(:,j,:), [round(axialShift_tilt_fit(j)), 0]);
end

disp('OCTA Processing complete');

end
