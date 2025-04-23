function [avgOCT_mcorr, OCTA_mcorr, yShift_axial] = mcorrGlobal(avgOCT_lin, OCTA, usfac, limit)

%%% Global motion parameter %%%
% limit = 30;
% usfac = 1; 

numFrames = size(avgOCT_lin, 3);
yShift_axial = zeros([numFrames 1]);

for I = 1:numFrames
    %%% Every 'for' loop, reference frame will be the middle frame %%%
    [output, ~] = dftregistration(fft2(avgOCT_lin(:, :, round(numFrames./2))), fft2(avgOCT_lin(:, :, I)), usfac);
    
    %%% Assign and save the shifting value for lateral (xShift) and axial (yShift) %%%
    yShift_axial(I) = round(output(3));
    
    %%% Thresholding  value was found via plotting the shifting values %%%
    if abs(output(3)) >= limit
        output(3) = 0;
    end
    
    avgOCT_mcorr(:,  :, I)  = circshift(avgOCT_lin(:, :, I), [output(3), 0]);
    OCTA_mcorr(:,  :, I) = circshift(OCTA(:, :, I), [output(3), 0]);
end
subplot(1,2,1),imagesc(imadjust(mat2gray(avgOCT_mcorr(:,:,round(numFrames./2))))),colormap(gray),title('OCT')
subplot(1,2,2),imagesc(imadjust(mat2gray(OCTA_mcorr(:,:,round(numFrames./2))))),colormap(gray),title('OCTA')

% OCT enhacement for segmentation %
avgOCT_enhance = avgOCT_mcorr;
for i = 1:size(avgOCT_mcorr,3)
    if i == 1
        avgOCT_enhance(:,:,i) = mean(avgOCT_mcorr(:,:,i:i+1),3);
    elseif i == size(avgOCT_mcorr,3)
        avgOCT_enhance(:,:,i) = mean(avgOCT_mcorr(:,:,i-1:i),3);
    else
        avgOCT_enhance(:,:,i) = mean(avgOCT_mcorr(:,:,i-1:i+1),3);
    end
end
