% Display DOPU image %
clearvars
addpath('C:\Users\KoreanBorg\Desktop\Reseachcore\Project-Melanoma\code');

disp('Filtering DOPU...');
DOPU_filt=DOPU;
DOPU_filt(DOPU_filt>0.95) = 1; % threshold the DOPU
DOPU_test=medfilt3(DOPU_filt,[3 5 3]); %[3 5(or 3) 3] (height, wdith, depth) with B-scans, [3 5 13] without B-scans
disp('DOPU filtered.');

% imagesc(DOPU_Bscans(:,:,1,500)),colormap('hot')
for i = 1:size(DOPU_test,3)
    imagesc(flipud(DOPU_test(:,:,i)),[0,1]);colormap(cmap_dopu_r)
    pause(0.1)
end

% DOPU B-scan
ii=250;
imgD_rgb_r = DOPU_Bscans(:,:,:,ii);
[imgD_ind, cmap]=rgb2ind(imgD_rgb_r,32);
figure;imagesc(imgD_ind);colormap(cmap);title('DOPU B-scan')
