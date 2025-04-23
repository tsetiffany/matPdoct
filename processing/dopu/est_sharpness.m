function [B, A, F] = est_sharpness(img)
% Estimate image sharpness usign three different methods

% BRISQUE model calculates image quality, lower equals less distortion
B = brisque(img);

% Measuring accutance (edge contrast of image) using the mean of the
% graident filter (direction doesn't matter, only mangitude)
[Gmag, ~] = imgradient(img);
A = mean(Gmag(:));

% Frequency analysis of the image; higher frequencies means a sharper image
% in Matalb, unless performing an fftshift, the highest frequency values
% are at the center of the FFT

%sum up the magnitude to find what amoutn of the total is within the top
%20% of the frequencies
[height, width] = size (img);
stepH = round(height/5);
stepW = round(width/5);
imgFFT = fft2(img);
imgFFT_highfreq = abs(imgFFT(stepH*2:stepH*3,stepW*2:stepW*3));
F = sum(imgFFT_highfreq(:))/sum(abs(imgFFT(:)));

end