%% DOPU B-scan overlay
% requires DOPU_test, avgOCT, cmap_dopu_r
% saves matrix of DOPU B-scans (in rgb), a tiff stack, and 
%% DOPU combined B-scans with log scale avgOCT
numPoints = size(DOPU,1);
numAlines = size(DOPU,2);
numBscans = size(DOPU,3);

% prepare matalb matrix for the Bscans in rgb
DOPU_Bscans = zeros(numPoints, numAlines, 3, numBscans); %matrix for the rgb images (index with the frame number LAST!)

% prepare Tiff object
% t = Tiff(fullfile(loadloc,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'.tif']),'w');
outputFileName = fullfile(loadloc,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'.tif']);

% loop over frames
for ii=1:size(DOPU,3)
    % ii= 200;
    
    % preparing the greyscale avgOCT luminance mask
    
    %for when there are zeros in the avgOCT (set them to around 13-30)
    G = flipud(squeeze(avgOCT(:,:,ii)));
    G(G==0)=30;
    imgL=10*log10(G); 
    
%     imgL=10*log10(flipud(avgOCT(:,:,ii))); % log scale tthe avgOCT (is already double, so only need 10log10 not 20log10)
    imgL=mat2gray(imgL);
%     imgL_test=imadjust(imgL,[0,1],[],2.2);
    imgL_test=imadjust(imgL,[0.4,0.8],[]); % perform some imadjustment
%     figure;subplot(1,2,1);imagesc(imgL);colormap(gray);
%     subplot(1,2,2);imagesc(imgL_test);colormap(gray);

    % convert the DOPU image to RGB then HSV
    imgD=flipud(DOPU_test(:,:,ii));
    colorMapForSaving = imresize(cmap_dopu_r,[256, 3], 'nearest');
    imgD = uint8(255.0 *mat2gray(imgD,[0,1])); %without mat2gray the images lose resolution and look like gray2ind, even though mat2gray on ingD doesn't seem to change teh values
    %     imgD = uint8(255.0 *imgD); 
    % imgD_rgb=ind2rgb(floor(imgD.*255),colorMapForSaving);
    imgD_rgb=ind2rgb(floor(imgD),colorMapForSaving);
    imgD_hsv=rgb2hsv(imgD_rgb); %h=1, s=2, v=3
    
    % replace the V channel with the log image
    imgD_hsv(:,:,3)=imgL_test;
    imgD_rgb_r=hsv2rgb(imgD_hsv); % final combined B-scan (is rgb so a 3D matrix)
    
%     %%%%%%%%%% COMMENT OUT WHEN RUNNING!
%     figure;
%     ax1=subplot(2,2,1);imshow(imgD);colormap(ax1,cmap_dopu_r);
%     subplot(2,2,2);imshow(imgD_rgb);
%     ax3=subplot(2,2,3);imagesc(imgL_test);colormap(ax3,gray);
%     subplot(2,2,4);imshow(imgD_rgb_r);
    
    
    % save into a matrix
    DOPU_Bscans(:,:,:,ii) = imgD_rgb_r;
    
    % save as a tiff stack
%     write(t,imgD_rgb_r);
    imwrite(imgD_rgb_r, outputFileName, 'WriteMode', 'append',  'Compression','none');
    
    fprintf('Frame processed : %d\n', ii);
end
% close(t);

figure;
ax1=subplot(2,2,1);imagesc(flipud(DOPU_test(:,:,ii)));colormap(ax1,cmap_dopu_r);
subplot(2,2,2);imshow(imgD_rgb);
ax3=subplot(2,2,3);imagesc(imgL_test);colormap(ax3,gray);
subplot(2,2,4);imshow(imgD_rgb_r);

save(fullfile(loadloc,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'.mat']), 'DOPU_Bscans', '-v7.3');

%convert to indexed image for display
figure;
for ii=1:size(DOPU,3)
    imgD_rgb_r = DOPU_Bscans(:,:,:,ii);
    [imgD_ind, cmap]=rgb2ind(imgD_rgb_r,56); % 64 is the highest number allowable if the bg is black
    imagesc(imgD_ind);colormap(cmap);
    pause(0.05);
end


% imgComb=(1-imgL_test).*imgD;
% figure;imagesc(imgComb);colormap(cmap_dopu);
