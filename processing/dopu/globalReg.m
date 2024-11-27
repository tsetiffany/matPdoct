function volume_mcorr = globalReg(volume, usfac, ref, mcorrpath, fpsave, volNo)


%% Introduction for Destiny %%
% The final output will be volume_mcorr and you can run who frame and take a look

volume = avgOCT; % Try just use any type of avgOCT or OCTA and change it into volume


numFrames = size(volume, 3); % Catch the max numframes value


usfac = 100; % This is highest compensation factor usually Francis put it 100

ref = volume(: , : , 20); % Take out any reference frame to match every frame of B-scan

% You can ignore this code %
% ProcessBar = waitbar( 0 , 'Preparing the process'); % Just I am using this to see how the process is going 
% You can ignore this code %

for ii = 1:numFrames
    [output, ~] = dftregistration(fft2(ref), fft2(imgaussfilt(abs(volume(:, :, ii)), 2)), usfac);
    yShift(ii) = round(output(3));
%    xShift(ii) = round(output(4));
    volume_mcorr(:, :, ii) = circshift(volume(:, :, ii), yShift(ii), 1);
%    volume_mcorr(:, :, ii) = circshift(volume(:, :, ii), [yShift(ii) xShift(ii)]);


%     BarValue = round(ii./numFrames*100);
    
%     waitbar(BarValue/100 , ProcessBar , sprintf('%d% has been completed.. Stand by...', BarValue));
    
end

% save(fullfile(mcorrpath,[fpsave,'_',num2str(volNo-1),'    _mcorr.mat']), 'volume_mcorr', '-v7.3');



%% To check before and after the globalReg %%

figure('Position' , [200 200 1500 600]) 

for i = 1:length(volume_mcorr)
    OCT = 20.*log10(abs(squeeze(OCT_LB(: , : , i))));
    OCT_New = 20.*log10(abs(squeeze(OCT_SB(: , : , i))));
    
    
    
%     OCT_Plot = abs(hilbert(OCT_LB(: , 234 , i)));
%     OCT_New_Plot = abs(hilbert(OCT_SB(: , 282 , i)));
    
    

    
subplot(1 , 2 , 1);
imagesc(imadjust(mat2gray(volume_mcorr)))); colormap(gray) ;
title('After the correction');
      
    
    
    % Put your original version of OCTA or avgOCT or volume_mcorr %
    subplot(1 , 2 , 2);  
    imagesc(imadjust(mat2gray(avgOCT))); colormap(gray); 
    title('Before the correction');
    
    
    pause(0.005)
end




end

