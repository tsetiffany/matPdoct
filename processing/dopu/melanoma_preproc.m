% Melanoma pre-processing from ProcdData data; separating the channels and
% performing 

close all; clear all; clc;

%% set up loading directories
% loadloc = 'C:\Users\Destiny\Documents\Research\MATLAB\Melanoma';    % directory containing the volume(s)
loadloc = 'F:\Clinical Data\2021.11.09\2021.11.09_P\12_04_48'; 
addpath(loadloc);
savepath = loadloc;                 % directory to save the files to
script_dir = 'C:\Users\KoreanBorg\Desktop\Reseachcore\Project-PD_OCT\Melanoma\code\JonesMatrixCode\JonesMatrixCode_Des_Ver6_ProcForPSSAO-OCT';
addpath(script_dir);                % directory containing this script

%%%% Set parameters %%%%
numBMscans = 1;     % number of BM scans
numFrames = 1000;   % total number of frames (including all BM scans)
strtFrame = 1;      % which frame to start processing (typically 1 for most purposes)

%% Load the file(s) %%%
filename = '12_04_48-' ; %%%%%%%%%%%%% Change depending on file name %%%%%%%%%%%%
disp(filename);
fn_num=filename;
% disp('Loading file...');
% load(fullfile(loadloc,filename));
% disp('File loaded');
% 
% % Extract sizes
% [~,numAlines,numBscans] = size(ProcdData);
% 
disp('Loading files...');
load(fullfile(loadloc,[filename,'A.mat']));
load(fullfile(loadloc,[filename,'B.mat']));

% Load colormap %
colormap_dir = 'C:\Users\KoreanBorg\Desktop\Reseachcore\Project-PD_OCT\Melanoma\code\JonesMatrixCode\data';
load(fullfile(colormap_dir,'cmap_dopu_r.mat'));    %load the DOPU colormap
load(fullfile(colormap_dir,'cmap_RPE_r.mat'));     % load the RPE colourmap for en face projection of low DOPU values
load(fullfile(colormap_dir,'cmap_OCTA_r2.mat'));   % load the depth-encoded OCTA map

% %% Volume flythrough to set cropping boundaries
% figure;
% for ii=1:numBscans
%     img = mat2gray(squeeze(20.*log10(abs(ProcdData(:,:,ii)))));
%     imagesc(imadjust(img)); colormap(gray);title(ii);
%     pause(0.02);
% end
ref_Frame = 500;
figure(1),imagesc(imadjust(mat2gray(20*log10(abs(cplxData_A(:,:,ref_Frame)))))),colormap('gray')
%% Reference frame process%%
refOCT_P = cplxData_A(:,:,ref_Frame);
refOCT_S = cplxData_B(:,:,ref_Frame);
refPhaseOff = repmat(angle(sum(refOCT_S.*conj(refOCT_P),1)), [size(refOCT_S,1) 1]); % new paper; get the value from is/os only? check later
refOCT_Cplx = refOCT_P + (refOCT_S.*exp(-1j.*refPhaseOff));

% Noise estimate %
for i = 1:size(refOCT_P,2)
    nOCT_P_real = var(real(refOCT_P(end-51:end,i))); % orig 1:20
    nOCT_P_imag = var(imag(refOCT_P(end-51:end,i)));
    nOCT_P_cplx(i) = nOCT_P_real + nOCT_P_imag; % Should this be actualy complex noise?????
    nOCT_S_real = var(real(refOCT_S(end-21:end,i)));
    nOCT_S_imag = var(imag(refOCT_S(end-21:end,i)));
    nOCT_S_cplx(i) = nOCT_S_real + nOCT_S_imag;
end
OCT_PN = median(nOCT_P_cplx)+std(nOCT_P_cplx);
OCT_SN = median(nOCT_S_cplx)+std(nOCT_P_cplx);

% Cropping
depthROI = [51 650];
for ii = 1:25:size(cplxData_A,3)
    imagesc(imadjust(mat2gray(20*log10(abs(cplxData_A(:,:,ii)))))),colormap('gray'), hold on
    plot(depthROI(1)*ones(1,size(cplxData_A,2)),'y','LineWidth',3),
    plot(depthROI(2)*ones(1,size(cplxData_A,2)),'y','LineWidth',3),
    pause(0.1)
end

%% Motion correction %%
cplxData_A_crop=cplxData_A(depthROI(1):depthROI(2),:,:);
cplxData_B_crop=cplxData_B(depthROI(1):depthROI(2),:,:);
ProcdData_ChP=cplxData_A_crop;
ProcdData_ChS=cplxData_B_crop;
clear cplxData_A_crop cplxData_B_crop
% ProcdData_ChP = ProcdData_ChP(ctop:cbot,:,:);
% ProcdData_ChS = ProcdData_ChS(ctop:cbot,:,:);

% reset sizes
[numPoints,numAlines,numBscans] = size(ProcdData_ChP);

% Motion correction
% addpath('/project/6007991/borg/STUDENTS/destinyh/MATLAB/volume_mcorr_code');
[~, yShift_global, ~] = MotionCorrectionGlobal(ProcdData_ChP, filename, savepath);

volume_mcorr_ChP = zeros(size(ProcdData_ChP));
volume_mcorr_ChS = zeros(size(ProcdData_ChS));

for I=1:numBscans
    volume_mcorr_ChS(:,  :, I) = circshift(ProcdData_ChS(:, :, I), [yShift_global(I), 0]);
    volume_mcorr_ChP(:,  :, I) = circshift(ProcdData_ChP(:, :, I), [yShift_global(I), 0]);
end

figure;plot(yShift_global);pause(1);

figure;
for ii=1:10:numBscans
    subplot(2,2,1)
    imgP = mat2gray(squeeze(20.*log10(abs(ProcdData_ChP(:,:,ii)))));
    imagesc(imadjust(imgP)); colormap(gray);title('ChP');xlabel(ii);
    
    subplot(2,2,2)
    imgPM = mat2gray(squeeze(20.*log10(abs(volume_mcorr_ChP(:,:,ii)))));
    imagesc(imadjust(imgPM)); colormap(gray);title('ChP - mcorr');xlabel(ii);
    
    subplot(2,2,3)
    imgS = mat2gray(squeeze(20.*log10(abs(ProcdData_ChS(:,:,ii)))));
    imagesc(imadjust(imgS)); colormap(gray);title('ChS');xlabel(ii);

    subplot(2,2,4)
    imgSM = mat2gray(squeeze(20.*log10(abs(volume_mcorr_ChS(:,:,ii)))));
    imagesc(imadjust(imgSM)); colormap(gray);title('ChS - mcorr');xlabel(ii);
    
    pause(0.01);
end


% Central slow scan
% ii=500;
% figure;
% subplot(2,2,1)
% imgP = mat2gray(squeeze(20.*log10(abs(ProcdData_ChP(:,ii,:)))));
% imagesc(imadjust(imgP)); colormap(gray);title('ChP');
% 
% subplot(2,2,2)
% imgPM = mat2gray(squeeze(20.*log10(abs(volume_mcorr_ChP(:,ii,:)))));
% imagesc(imadjust(imgPM)); colormap(gray);title('ChP - mcorr');
% 
% subplot(2,2,3)
% imgS = mat2gray(squeeze(20.*log10(abs(ProcdData_ChS(:,ii,:)))));
% imagesc(imadjust(imgS)); colormap(gray);title('ChS');
% 
% subplot(2,2,4)
% imgSM = mat2gray(squeeze(20.*log10(abs(volume_mcorr_ChS(:,ii,:)))));
% imagesc(imadjust(imgSM)); colormap(gray);title('ChS - mcorr');

%% More cropping if necessary
% Cropping
depthROI = [51 600];
for ii = 1:25:size(volume_mcorr_ChP,3)
    imagesc(imadjust(mat2gray(20*log10(abs(volume_mcorr_ChP(:,:,ii)))))),colormap('gray'), hold on
    plot(depthROI(1)*ones(1,size(volume_mcorr_ChP,2)),'y','LineWidth',3),
    plot(depthROI(2)*ones(1,size(volume_mcorr_ChP,2)),'y','LineWidth',3),
    pause(0.1)
end

%% saving mcorr files
volume_mcorr_ChP = volume_mcorr_ChP(depthROI(1):depthROI(2),:,:);
volume_mcorr_ChS = volume_mcorr_ChS(depthROI(1):depthROI(2),:,:);

noise_PS = [OCT_PN, OCT_SN];

disp('Saving mcorrs...');
save(fullfile(savepath,[filename,'    _mcorr_ChP_Global.mat']), 'volume_mcorr_ChP', '-v7.3');
save(fullfile(savepath,[filename,'    _mcorr_ChS_Global.mat']), 'volume_mcorr_ChS', '-v7.3');
disp('mcorrs saved');
disp('Saving noise...');
save(fullfile(savepath,[filename,'    _noise.mat']), 'noise_PS', '-v7.3');
disp('noise saved');

%% DOPU volume processing %%
tic;
iF=1; % frame counter
OCT_PN = noise_PS(1);
OCT_SN = noise_PS(2);

%%% Processing options %%%
adapt = 0;          % 1 (adaptive averaging kernel) / 0 (rigid kernel)
raster = 0;         % 1 (raster scan, requires some flyback cropping) / 0 (bi-directional)
kernel = [3,5];     % averaging kernel dimensions (axial, lateral)

clearvars dopu CompCplx
for I = 1:numBMscans:numFrames
    
    K = ((I-1)/numBMscans)+1;
    % JM % 
    for J = 1:numBMscans
        OCT_P  = volume_mcorr_ChP(:,:,(I+J-1));  
        OCT_S  = volume_mcorr_ChS(:,:,(I+J-1));
        S0    = OCT_P.*conj(OCT_P) + OCT_S.*conj(OCT_S);
        S1    = OCT_P.*conj(OCT_P) - OCT_S.*conj(OCT_S);
        S2    = 2.*real(OCT_P.*conj(OCT_S));
        S3    = 2.*imag(OCT_P.*conj(OCT_S));
        S0_NC = S0 - (OCT_PN + OCT_SN);
        S1_NC = S1 - (OCT_PN - OCT_SN);

        if adapt == 1
            mS0   =  adaptKernelSmoothing_r(S0_NC, kernel,5,0);
            mS1   =  adaptKernelSmoothing_r(S1_NC, kernel,5,0);
            mS2   =  adaptKernelSmoothing_r(S2, kernel,5,0);
            mS3   =  adaptKernelSmoothing_r(S3, kernel,5,0);
        else
            mS0   =  smooth2DFilter(S0_NC, kernel);
            mS1   =  smooth2DFilter(S1_NC, kernel);
            mS2   =  smooth2DFilter(S2, kernel);
            mS3   =  smooth2DFilter(S3, kernel);
        end

        dopu_Numer  = sqrt(mS1.^2 + mS2.^2 + mS3.^2);
        dopu_Denom  = mS0;
        dopu(:,:,J) = dopu_Numer./dopu_Denom; 
        rPhaseOff = repmat(angle(sum(OCT_S.*conj(OCT_P),1)), [size(OCT_S,1) 1]);
        OCT_Cplx = (OCT_P) + (OCT_S.*exp(-1j.*rPhaseOff));
%             Xconj    = OCT_Cplx.*conj(refOCT_Cplx);
%             bPhaseOff  = repmat(angle(sum(Xconj,1)), [size(Xconj,1) 1]);
%             CompCplx(:,:,J) = OCT_Cplx.*exp(-1j*bPhaseOff);
        CompCplx(:,:,J) = OCT_Cplx;
    end

%         figure(1),
%         subplot(1,3,1),imagesc(imadjust(mat2gray(20*log10(abs(CompCplx))))),colormap('gray')
%         subplot(1,3,2),imagesc(imadjust(mat2gray(((dopu))))),colormap('gray')
%         DOPU_filt=abs(dopu);
%         DOPU_filt(DOPU_filt>0.95) = 1; % threshold the DOPU
%         DOPU_test=medfilt3(DOPU_filt,[3 5 3]);
%         subplot(1,3,3),imagesc((DOPU_test),[0,1]);colormap(cmap_dopu_r);

%         mDOPU          = abs(10.*log10(dopu)); % used this for the 1 BM
%         scan mouse data, because otherwise the DOPU looked like garbage

    if numBMscans == 1
        avgOCT(:,:,K)  = abs(CompCplx);
        mDOPU          = dopu;
        OCTA = [];
    else
        avgOCT(:,:,K)  = mean(abs(CompCplx),3);
        mDOPU          = mean(dopu, 3);
        if numBMscans == 2
            OCTA(:,:,K) = abs(CompCplx(:,:,2)-CompCplx(:,:,1));
        else
            OCTA(:,:,K) = var(CompCplx, 0, 3);
        end
    end

    mDOPU          = dopu;
    mDOPU          = 1./mDOPU;
    mDOPU(mDOPU<1) = 1;
    DOPU(:,:,K)   = 1./mDOPU;

    DOPU_filt=1./mDOPU;
    DOPU_filt(DOPU_filt>0.95) = 1; % threshold the DOPU
    DOPU_test=medfilt3(DOPU_filt,[3 5 3]);
    imagesc((DOPU_test),[0,1]);colormap(cmap_dopu_r);
    
    fprintf('Frame processed : %d\n', iF);
    iF = iF+1;
end
toc
save(fullfile(savepath,[fn_num,'_avgOCT.mat']), 'avgOCT', '-v7.3');

% reset sizes
[numPoints,numAlines,numBscans] = size(DOPU);
%% Display the DOPU images %%
disp('Filtering DOPU...');
DOPU_filt=DOPU;
DOPU_filt(DOPU_filt>0.95) = 1; % threshold the DOPU
DOPU_test=medfilt3(DOPU_filt,[3 5 3]); %[3 5(or 3) 3] (height, wdith, depth) with B-scans, [3 5 13] without B-scans
disp('DOPU filtered.');

% ii=320;
for ii = 301:700
    imgD=(DOPU_test(:,:,ii));
    figure(10);
    imagesc(imgD,[0,1]);colormap(cmap_dopu_r);colorbar;
    if adapt==1
        title(['DOPU, adaptive, [',num2str(kernel(1)),',',num2str(kernel(2)),']']);
    else
        title(['DOPU, rigid, [',num2str(kernel(1)),',',num2str(kernel(2)),']']);
    end
    xlabel(ii);
end
if adapt == 1
    save(fullfile(savepath,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'_adaptSmooth_dopu.mat']), 'imgD', '-v7.3');
else
    save(fullfile(savepath,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'_regid_dopu.mat']), 'imgD', '-v7.3');
end

%% DOPU combined B-scans with log scale avgOCT
% prepare matalb matrix for the Bscans in rgb
DOPU_Bscans = zeros(numPoints, numAlines, 3, numBscans); %matrix for the rgb images (index with the frame number LAST!)

% prepare Tiff object
outputFileName = fullfile(loadloc,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'.tif']);

% loop over frames
for ii=1:size(DOPU,3)
    
    % preparing the greyscale avgOCT luminance mask
    if ii==1 || ii==size(DOPU,3)
        G = flipud(squeeze(avgOCT(:,:,ii)));
    else
        % average one frame prior and after
        G = flipud(mean(avgOCT(:,:,ii-1:ii+1),3));
    end
    %for when there are zeros in the avgOCT (set them to around 13-30)
    G(G==0)=30;
    imgL=10*log10(G); 
    
    imgL=mat2gray(imgL);
    imgL_test=imadjust(imgL,[0.25,0.8],[]); % perform some imadjustment
    % ******************* COMMENT OUT BEFORE RUNNING!!! **************
%     figure;subplot(1,2,1);imagesc(imgL);colormap(gray);
%     subplot(1,2,2);imagesc(imgL_test);colormap(gray);

    % convert the DOPU image to RGB then HSV
    imgD=flipud(DOPU_test(:,:,ii));
    colorMapForSaving = imresize(cmap_dopu_r,[256, 3], 'nearest'); % interpolate to 256 colour values (one for eevry number possible in the image)
    imgD = uint8(255.0 *mat2gray(imgD,[0,1])); %without mat2gray the images lose resolution and look like gray2ind, even though mat2gray on ingD doesn't seem to change teh values
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
%     imwrite(imgD_rgb_r, outputFileName, 'WriteMode', 'append',  'Compression','none');
    
    fprintf('Frame processed : %d\n', ii);
end

if adapt == 1
    save(fullfile(savepath,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'_adaptSmooth_dopuBscan.mat']), 'DOPU_Bscans', '-v7.3');
else
    save(fullfile(savepath,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'_regid_dopuBscan.mat']), 'DOPU_Bscans', '-v7.3');
end

%% create movie file
% create the video writer with 1 fps
knum=sprintf('%d%02d',kernel(1),kernel(2));
if adapt==1
	vidname=fullfile(savepath,[fn_num,'adaptive',knum,'_montage.avi']);
else
    vidname=fullfile(savepath,[fn_num,'rigid',knum,'_montage.avi']);
end
writerObj = VideoWriter(vidname);
writerObj.FrameRate = 10; % set the seconds per image

% open the video writer
open(writerObj);

for ii = 1:size(avgOCT,3)-1
    figure(1)  
    
    % logscale avgOCT
     if ii==1 || ii==size(avgOCT,3)
        G = flipud(squeeze(avgOCT(:,:,ii)));
    else
        % average one frame prior and after for display
        G = flipud(mean(avgOCT(:,:,ii-1:ii+1),3));
     end
    G(G==0)=30;
    imgL=10*log10(G); 
    imgL=mat2gray(imgL);
    imgL_test=imadjust(imgL,[0.25,0.9],[]);
    
    % OCTA
%     imgA = mat2gray(flipud(OCTA(:,:,ii)));
%     imgA = imadjust(imgA,[0,0.2],[],2.2);
    imgA=[];
    
    % DOPU
    imgD=flipud(DOPU_test(:,:,ii));
    
    % DOPU B-scans
    imgD_rgb_r = DOPU_Bscans(:,:,:,ii);
    [imgD_ind, cmap]=rgb2ind(imgD_rgb_r,32); % 64 is the highest number allowable if the bg is black
    
    ax1=subplot(3,1,1);imagesc(imgL_test);colormap(ax1,gray);title('log avgOCT');xlabel(ii);colorbar
%     ax2=subplot(1,3,2);imagesc(imgA);colormap(ax2,gray);title('OCTA');xlabel(ii);
    ax2=subplot(3,1,2);imagesc(imgD_ind);colormap(ax2,cmap);title('DOPU B-scan');xlabel(ii);colorbar
    ax3=subplot(3,1,3);imagesc(imgD,[0,1]);colormap(ax3,cmap_dopu_r);colorbar;title('DOPU');xlabel(ii);
    set(gcf,'Position',[200 100 400 700])

    % convert the image to a frame
    frame = getframe(gcf) ;

    % write the frames to the video
    writeVideo(writerObj, frame);
    pause(0.01)
end

% close the writer object
close(writerObj);

%%%% (optional fundus)
f=(squeeze(mean(avgOCT,1)));
    f(f==0)=30;
    imgF=10*log10(f); 
    imgF=mat2gray(imgF);
    imgF_test=imadjust(imgF,[0.25,0.9],[]);
figure;imagesc(imgF);colormap(gray);axis equal; axis tight;
imwrite(imgF,fullfile(savepath,[filename,'_FUND.tif']));