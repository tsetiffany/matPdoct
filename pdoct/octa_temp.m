%% OCTA Processing Pipeline %%
%% Single Channel PDOCT 
%% Update: 06 December 2024 (MSI)

% Load the volume and check flythrough
figure(),
for i = 1:size(cplxData_A,3)
    imagesc(imadjust(mat2gray(20*log10(abs(cplxData_A(:,:,i))))));colormap(gray),title(sprintf('%d',i))
    pause(0.001)
end


%% Delete every first BScan from the BMscan set
indices_to_keep = setdiff(1:size(CompCplx, 3), 1:4:size(CompCplx, 3)); % 4BM
cplxData_A_1 = CompCplx(:, :, indices_to_keep);
cplxData_A_1 = cplxData_A_1(:, 15:485,:); % cut flyback

%Location of motion correction code
addpath('C:\Users\tiffa\Documents\1. Projects\PDOCT Processing\main_vctrl\Code_workingversion\matPdoct\MotionCorrection');

%%% Local sub-pixel motion correction %%%
numMscans = 3;
S_Mcorr  = cplxData_A_1;
numFrames = size(cplxData_A_1, 3);
for I = 1:numMscans:numFrames-1
    batchFrames = S_Mcorr(:, :, I:I+numMscans-1);
    correctedBatch = MotionCorrection(batchFrames);
    S_Mcorr(:, :, I:I+numMscans-1) = correctedBatch;
end

clearvars -except S_Mcorr 


%% OCTA
volume = S_Mcorr;

numMscans = 3;
for I = 1:numMscans:size(volume,3)
		K = ((I-1)/numMscans)+1;
		
		Xconj_2   = volume(:,:,I+1).*conj(volume(:,:,I));
		Xconj_3   = volume(:,:,I+2).*conj(volume(:,:,I));
%         Xconj_4   = volume(:,:,I+3).*conj(volume(:,:,I));

		BulkOff_2 = repmat(angle(sum(Xconj_2)), [size(Xconj_2,1) 1]);
		BulkOff_3 = repmat(angle(sum(Xconj_3)), [size(Xconj_3,1) 1]);
%         BulkOff_4 = repmat(angle(sum(Xconj_4)), [size(Xconj_4,1) 1]);

		Bscan_1  = volume(:,:,I);
		Bscan_2  = volume(:,:,I+1) .* exp(-1j*BulkOff_2);
		Bscan_3  = volume(:,:,I+2) .* exp(-1j*BulkOff_3);
%         Bscan_4  = volume(:,:,I+3) .* exp(-1j*BulkOff_3);

        % Average OCT %
% 		avgOCT(:,:,K) = (abs(Bscan_1) + abs(Bscan_2) + abs(Bscan_3) + abs(Bscan_4))/4;
        avgOCT(:,:,K) = (abs(Bscan_1) + abs(Bscan_2) + abs(Bscan_3))/3;


		% Variance %
% 		Var(:,:,K) = abs(var(cat(3,(Bscan_1),(Bscan_2),(Bscan_3)),0,3));
%       Sub(:,:,K) = abs(Bscan_1 - Bscan_2) + abs(Bscan_2 - Bscan_3) + abs(Bscan_3 - Bscan_4);
        Sub(:,:,K) = abs(Bscan_1 - Bscan_2) + abs(Bscan_2 - Bscan_3);

		fprintf('OCTA volume process: %d\n', K);
end

%Check
figure(),
for i = 1:size(Var,3)
    imagesc(imadjust(mat2gray(abs(Var(:,:,i)))));colormap(gray),title(sprintf('%d',i))
    pause(0.001)
end

clearvars -except avgOCT Sub


%% Motion Correction

% Global Motion
numFrames   = size(avgOCT,3);
yShift_axial = zeros([numFrames 1]);
limit = 150;
usfac = 1; 
avgOCT_mcorr = avgOCT;
OCTA_mcorr_Var = Sub;
% OCTA_mcorr_Var = Var_interleaved;

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

end

% Global motion correction based on smoothen avgOCT motion correction
for I = 1:numFrames

    avgOCT_mcorr(:,:,I)  = circshift(avgOCT(:,:,I), [yShift_axial(I), 0]);
    OCTA_mcorr_Var(:,:,I)  = circshift(Sub(:,:,I), [yShift_axial(I), 0]);
end


% Tilt correction
numLines = size(avgOCT_mcorr,2);
usfac = 1;

avgOCT_tcorr        = zeros(size(avgOCT_mcorr));
OCTA_tcorr_Var      = zeros(size(OCTA_mcorr_Var));

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

    %Motion correction processing as shifting value along axial direction
    avgOCT_tcorr(:,j,:) = circshift(avgOCT_mcorr(:,j,:), [round(axialShift_tilt_fit(j)), 0]);
    OCTA_tcorr_Var(:,j,:) = circshift(OCTA_mcorr_Var(:,j,:), [round(axialShift_tilt_fit(j)), 0]);
    % OCTA_tcorr_var(:,j,:) = circshift(var_mcorr(:,j,:), [round(axialShift_tilt_fit(j)), 0]);

end

OCTA_tcorr_sub_1 = flip(OCTA_tcorr_Var);




%% Checking the volume

figure(),
for i = 1:size(OCTA_tcorr_Var,3)
    imagesc(imadjust(mat2gray(abs(OCTA_tcorr_Var(:,:,i)))));colormap(gray),title(sprintf('%d',i))
    pause()
end

figure(),
for i = 250:size(OCTA_tcorr_Var,1)
    imagesc(imadjust(mat2gray(squeeze(OCTA_tcorr_Var(i,:,:)))));colormap(gray),title(sprintf('%d',i))
    pause()
end

roi = 320:327;
enf = mat2gray(squeeze(max(OCTA_tcorr_Var(roi,:,:),[],1)));

imagesc(imadjust(enf)); colormap gray;


% figure(),
% for i = 750:size(OCTA_tcorr_sub_1,1)
%     imagesc(imadjust(mat2gray(squeeze(mean(abs(OCTA_tcorr_sub_1(i,:,:)),1)))));colormap('gray'),title(sprintf('%d',i))
%     pause()
% end


figure(),imagesc(imadjust(mat2gray(squeeze(max(abs(OCTA_tcorr_sub_1(570:640,:,:)),[],1)))));colormap('gray')
% 
% 
% 
% ROI = 780:785;
% 
% enf_ini = mat2gray(squeeze(max(abs(OCTA_tcorr_sub_1(ROI,:,:)),[],1)));
% enf_hist2 = adapthisteq(enf_ini,"NumTiles",[25 25],"ClipLimit",0.002,"Distribution","uniform");
% enf_enh2 = imcenhance(enf_hist2);
% enf_filt2 = medfilt2(enf_enh2,[3 3]);
% figure(),imagesc(imadjust(enf_filt2)); colormap gray