function [volume_mcorr_PS, volume_mcorr_ChP, volume_mcorr_ChS] = globalReg_PS(volume_PS, volume_ChP, volume_ChS, usfac , ref)

%% Introduction for Destiny %%
% The final output will be volume_mcorr and you can run who frame and take a look

% volume = volume_mcorr; % Try just use any type of avgOCT or OCTA and change it into volume
numFrames = size(volume_PS, 3); % Catch the max numframes value
% usfac = 100; % This is highest compensation factor usually Francis put it 100
% ref = volume(: , : , 20); % Take out any reference frame to match every frame of B-scan


% You can ignore this code %
% ProcessBar = waitbar( 0 , 'Preparing the process'); % Just I am using this to see how the process is going 
% You can ignore this code %

for ii = 1:numFrames
    
    [output, ~] = dftregistration(fft2(ref), fft2(imgaussfilt(abs(volume_PS(:, :, ii)), 2)), usfac);
    yShift(ii) = round(output(3));
%    xShift(ii) = round(output(4));
    volume_mcorr_PS(:, :, ii) = circshift(volume_PS(:, :, ii), yShift(ii), 1);
    volume_mcorr_ChP(:, :, ii) = circshift(volume_ChP(:, :, ii), yShift(ii), 1);
    volume_mcorr_ChS(:, :, ii) = circshift(volume_ChS(:, :, ii), yShift(ii), 1);

%    volume_mcorr(:, :, ii) = circshift(volume(:, :, ii), [yShift(ii) xShift(ii)]);


end

end
% 
