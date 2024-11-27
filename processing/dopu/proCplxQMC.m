close all; clear all; clc;

loadloc = 'F:\Clinical Data\2021.11.09\2021.11.09_P\12_04_48';    % directory containing the volume(s)

savepath = loadloc;                 % directory to save the files to

% savepath = '/project/6007991/borg/STUDENTS/destinyh/DATA/Melanoma';                 % directory to save the files to
% script_dir = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\JonesMatrixCode_Des_Ver6_ProcForPSSAO-OCT';
script_dir = 'C:\Users\KoreanBorg\Desktop\Reseachcore\Project-PD_OCT\Melanoma\code\JonesMatrixCode\JonesMatrixCode_Des_Ver6_ProcForPSSAO-OCT';
addpath(script_dir);                % directory containing this script
%%
% loadloc_r = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data';
loadloc_r = 'C:\Users\KoreanBorg\Desktop\Reseachcore\Project-PD_OCT\Melanoma\code\JonesMatrixCode\data';
load(fullfile(loadloc_r,'cmap_dopu_r.mat'));    %load the DOPU colormap
load(fullfile(loadloc_r,'cmap_RPE_r.mat'));     % load the RPE colourmap for en face projection of low DOPU values
load(fullfile(loadloc_r,'cmap_OCTA_r2.mat'));   % load the depth-encoded OCTA map


%%%% Set parameters %%%%
numBMscans = 1;     % number of BM scans
numFrames = 1000;   % total number of frames (including all BM scans)
strtFrame = 1;      % which frame to start processing (typically 1 for most purposes)
kernel = [3,5];     % averaging kernel dimensions (axial, lateral)

%%% Processing options %%%
PS = 1;             % 1 (PS mode) / 0 (Non-PS mode)
adapt = 0;          % 1 (adaptive averaging kernel) / 0 (rigid kernel)
raster = 0;         % 1 (raster scan, requires some flyback cropping) / 0 (bi-directional)

%% Skip most of this if following melanoma_preproc directly!
%%% Load the files %%%
fn_num = '12_04_48-' ; %%%%%%%%%%%%% Change depending on file name %%%%%%%%%%%%
filename = fn_num;
disp(fn_num);

if PS == 1  
    fn_ChP = [fn_num, '    _mcorr_ChP_Global'];
    load(fullfile(loadloc,fn_ChP));
    if raster == 1
        cplxVol_ChP = volume_mcorr_ChP(:,51:end-50,:);
    else 
        cplxVol_ChP = volume_mcorr_ChP; %%%%%%%%%%%%% Only run this line after preproc
    end
    clear volume_mcorr_ChP %%%%%%%%%%%%%%%%%%%%5 Only run this line after preproc
    disp('P-channel mcorr loaded');
    
    fn_ChS = [fn_num,'    _mcorr_ChS_Global'];
    load(fullfile(loadloc,fn_ChS));
    if raster == 1
        cplxVol_ChS = volume_mcorr_ChS(:,51:end-50,:);
    else 
        cplxVol_ChS = volume_mcorr_ChS;%%%%%%%%%%%%% Only run this line after preproc
    end
    clear volume_mcorr_ChS%%%%%%%%%%%%% Only run this line after preproc
    disp('S-channel mcorr loaded');
else
    fn = [fn_num,'    _mcorr'];
    load(fullfile(loadloc,fn));
    if raster == 1
        cplxVol = volume_mcorr(:,51:end-50,:);
    else 
        cplxVol  = volume_mcorr;
    end
    clear volume_mcorr
    
    disp('mcorr loaded');
end
%% Size setup
% cplxVol_ChP=cplxData_A;
% cplxVol_ChS=cplxData_B;
% clear cplxData_A cplxData_B
% 
%%% size extraction and preallocating output spaces %%%
if PS == 1
    [numPoints,numAlines, ~] = size(cplxVol_ChS);
    DOPU = zeros(numPoints,numAlines,numFrames./numBMscans);
else
    [numPoints,numAlines, ~] = size(cplxVol);
end

numBscans = round(numFrames./numBMscans);

avgOCT = zeros(numPoints,numAlines,numBscans);
OCTA = zeros(numPoints,numAlines,numBscans);

% volumes are expected to be upside-down. Flip if they started rightside-up
cplxVol_ChS = flipud(cplxVol_ChS);
cplxVol_ChP = flipud(cplxVol_ChP);

for i =1:size(cplxVol_ChS,3)
    subplot(1,2,1),imagesc(imadjust(mat2gray(20*log10(abs(cplxVol_ChS(:,:,i))))))
    subplot(1,2,2),imagesc(imadjust(mat2gray(20*log10(abs(cplxVol_ChP(:,:,i)))))),colormap('gray')
    pause(0.1)
end

%% processing the files
tic;
iF=1; % frame counter
for I = 600%301:700%:numBMscans:numFrames
    K = ((I-1)/numBMscans)+1;
    if PS == 1 
        refOCT_P = cplxVol_ChP(:,:,300);
        refOCT_S = cplxVol_ChS(:,:,300);
        refPhaseOff = repmat(angle(sum(refOCT_S.*conj(refOCT_P),1)), [size(refOCT_S,1) 1]); % new paper; get the value from is/os only? check later
        refOCT_Cplx = refOCT_P + (refOCT_S.*exp(-1j.*refPhaseOff));
        
        % Noise estimate %
        for i = 1:size(refOCT_P,2)
            nOCT_P_real = var(real(refOCT_P(end-21:end,i))); % orig 1:20
            nOCT_P_imag = var(imag(refOCT_P(end-21:end,i)));
            nOCT_P_cplx(i) = nOCT_P_real + nOCT_P_imag; % Should this be actualy complex noise?????
            nOCT_S_real = var(real(refOCT_S(end-21:end,i)));
            nOCT_S_imag = var(imag(refOCT_S(end-21:end,i)));
            nOCT_S_cplx(i) = nOCT_S_real + nOCT_S_imag;
        end
        OCT_PN = median(nOCT_P_cplx);
        OCT_SN = median(nOCT_S_cplx);
        
        % JM % 
        for J = 1:numBMscans
            OCT_P  = cplxVol_ChP(:,:,(I+J-1));  
            OCT_S  = cplxVol_ChS(:,:,(I+J-1));
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
            Xconj    = OCT_Cplx.*conj(refOCT_Cplx);
            bPhaseOff  = repmat(angle(sum(Xconj,1)), [size(Xconj,1) 1]);
            CompCplx(:,:,J) = OCT_Cplx.*exp(-1j*bPhaseOff);
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
        
    else
        refOCT_Cplx = cplxVol(:,:,I);
        
        for J = 1:numBMscans
            OCT_Cplx    = cplxVol(:,:,(I+J-1));
            Xconj       = OCT_Cplx.*conj(refOCT_Cplx);
            bPhaseOff   = repmat(angle(sum(Xconj,1)), [size(Xconj,1) 1]);
            CompCplx(:,:,J) = OCT_Cplx.*exp(-1j*bPhaseOff);
        end
        
        avgOCT(:,:,K)  = mean(abs(CompCplx),3);
        if numBMscans == 2
            OCTA(:,:,K) = abs(CompCplx(:,:,2)-CompCplx(:,:,1));
        else
            OCTA(:,:,K) = var(CompCplx, 0, 3);
        end
    end
    
    fprintf('Frame processed : %d\n', iF);
    iF = iF+1;
end
toc


%% Save the images
disp('Saving files...');
if PS == 1
    if adapt == 1
        save(fullfile(savepath,[fn_num,'_DOPU_',num2str(kernel(1)),'_',num2str(kernel(2)),'_adaptSmooth_r.mat']), 'DOPU', '-v7.3');
    else
        save(fullfile(savepath,[fn_num,'_DOPU_',num2str(kernel(1)),'_',num2str(kernel(2)),'.mat']), 'DOPU', '-v7.3');
    end
end
save(fullfile(savepath,[fn_num,'_avgOCT.mat']), 'avgOCT', '-v7.3');
save(fullfile(savepath,[fn_num,'_OCTA.mat']), 'OCTA', '-v7.3');
disp('Files saved');


%%%%%%%%%%%%%  Display the images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Filtering DOPU...');
DOPU_filt=DOPU;
DOPU_filt(DOPU_filt>0.95) = 1; % threshold the DOPU
DOPU_test=medfilt3(DOPU_filt,[3 5 3]); %[3 5(or 3) 3] (height, wdith, depth) with B-scans, [3 5 13] without B-scans
disp('DOPU filtered.');


ii=320;
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
    save(fullfile(savepath,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'_adaptSmooth_r.mat']), 'DOPU_Bscans', '-v7.3');
else
    save(fullfile(savepath,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'.mat']), 'DOPU_Bscans', '-v7.3');
end

%% avgOCT, OCTA, DOPU and combined video montage

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
    
    ax1=subplot(2,2,1);imagesc(imgL_test);colormap(ax1,gray);title('log avgOCT');xlabel(ii);
    ax2=subplot(2,2,2);imagesc(imgA);colormap(ax2,gray);title('OCTA');xlabel(ii);
    ax3=subplot(2,2,3);imagesc(imgD_ind);colormap(ax3,cmap);title('DOPU B-scan');xlabel(ii);
    ax4=subplot(2,2,4);imagesc(imgD,[0,1]);colormap(ax4,cmap_dopu_r);colorbar;title('DOPU');xlabel(ii);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pigment and flow image (PAF aka QMC)
% load(fullfile(savepath,[strcat(filename,'_DOPU_3_5.mat')]));
% load(fullfile(savepath,[strcat(filename,'_avgOCT.mat')]));
% load(fullfile(savepath,[strcat(filename,'_OCTA.mat')]));


%% preparing OCTA & DOPU
% current process just imadjusts each B-scan, the previosu noise removal
% code was too heavy (added too much noise).  Should add histogram
% equalization across the frames as well.

% preparing OCTA
disp('Filtering OCTA...');
cplxOCTA=flipud(OCTA);
for ii=1:size(cplxOCTA,3)
    imgOCTA = mat2gray(cplxOCTA(:,:,ii));
    cplxOCTA(:,:,ii) = imadjust(imgOCTA,[0,0.2],[],2.2);
end
disp('OCTA filtered.');

% preparing DOPU
disp('Filtering DOPU...');
DOPU_filt=flipud(DOPU);

DOPU_filt(DOPU_filt>0.95) = 1; % threshold theDOPU
DOPU_test=(medfilt3(DOPU_filt, [3 5 3])); % [depth width frames]
disp('DOPU filtered.');

DOPU_test=flipud(DOPU_test);
% prepare the martices for the segmentation of the DOPU and set to zeros)
DOPU_depth_C = zeros(size(DOPU,2),size(DOPU,3));
DOPU_depth_C_bot = zeros(size(DOPU,2),size(DOPU,3));

%% Fundus image
disp('Generating fundus image...');

FUND = (imcenhance(mat2gray(squeeze(mean(avgOCT,1)))));

figure; imagesc(FUND);colormap gray; axis equal; axis tight; title('FUND');

imwrite(FUND,fullfile(savepath,[filename,'_FUND.tif']));

%% check DOPU image
botlim=5;
lowlim =10;
uplim = 20;

figure('pos',[50 50 1500 500]);
for ii=1:size(DOPU_test,3)-1 %cplxOCTA, DOPU_test
    ax1=subplot(1,2,1);
    imgA = squeeze(DOPU_test(:,:,ii));
    imagesc(imgA,[0,1]);colormap(ax1,cmap_dopu_r);colorbar;%colormap(gray, cmap_dopu);
    title(ii);
    hold on;
    plot(DOPU_depth_C(:,ii),'black','LineWidth', 2);
    plot(DOPU_depth_C_bot(:,ii),'black','LineWidth', 2);
    hold off;

%     ax2=subplot(1,2,2);
%     img2=mat2gray(flipud(squeeze(OCTA(:,:,ii))));
%     imagesc(img2, [0, 0.1]);colormap(ax2,gray);
%     title(ii);
%     
    hold on;
    plot(DOPU_depth_C(:,ii),'red','LineWidth', 1);
    plot(DOPU_depth_C(:,ii)-botlim,'blue','LineWidth', 1);
    plot(DOPU_depth_C(:,ii)-lowlim,'yellow','LineWidth', 1);
    plot(DOPU_depth_C_bot(:,ii),'green','LineWidth', 1);
    plot(DOPU_depth_C(:,ii)-uplim,'cyan','LineWidth', 1);
    
    hold off;
    
    pause(0.0001);
%     colorbar;
end


%% Finding the cut lines (segmentation lines) of the low DOPU region
disp('Preparing method C...');

% sizes shoudl already be set up earlier in code
% numPoints=size(DOPU_test,1);
% numAlines=size(DOPU_test,2);
% numBscans=size(DOPU_test,3);

startZ = 100; %depth to start searching at (usd to cut out noise)
thresh = 0.99; %threshold for finding low DOPU values (typically 0.95, lower values closer to LDR)
smoothing = 5;% if noisy, do 30, else do 5-20
botcrop = 10; % crop out any FPN near bottom

for B=1:numBscans % loop across the frames
    for A=1:numAlines-1 % loop across the A-lines
        
        % extract A-line (all depth values)
        Aline=squeeze(DOPU_test(:,A,B));
        
        % smooth out any outliers
        Aline=smoothdata(Aline,'gaussian',smoothing); 
        
        % find all indices where DOPU < thresh
        idx = find(Aline(startZ:end-botcrop)<thresh)+startZ;
        
        % extract the top and bottom values (if any are found), else set
        % the DOPU depth at the bottom of the volume
        if idx
            DOPU_depth_C(A,B)=idx(1);
            DOPU_depth_C_bot(A,B)=idx(end);
        else
            DOPU_depth_C(A,B)=numPoints;
            DOPU_depth_C_bot(A,B)=numPoints;
        end


    end
    % copy the second-last A-line just to avoid zero values
    DOPU_depth_C(end,:)=DOPU_depth_C(end-1,:);
    DOPU_depth_C_bot(end,:)=DOPU_depth_C_bot(end-1,:);
    
    % smooth out the final segmentation lines
    DOPU_depth_C(:,B)=smoothdata(DOPU_depth_C(:,B),'gaussian',50);
    DOPU_depth_C_bot(:,B)=smoothdata(DOPU_depth_C_bot(:,B),'gaussian',50);
    
    fprintf('DOPU frame: %d\n', B);
end

disp('Completed');

save(fullfile(savepath,[filename,'_DOPU_depth_C.mat']), 'DOPU_depth_C', '-v7.3');
save(fullfile(savepath,[filename,'_DOPU_depth_C_bot.mat']), 'DOPU_depth_C_bot', '-v7.3');

% displaying an elevation map
imgDDC = notchfilter(imcenhance(mat2gray((DOPU_depth_C))));%imcenhance (can do log10 on the dopu depth?

% display with inverted values(smaller numbers are higher elevation originally)
figure;i1=imagesc(1-imgDDC,[0,1]);title('DOPU_depth, method C','Interpreter', 'none');
colormap(hot);colorbar; axis equal; axis tight;

%%%%%%%%%%% check via the 'check DOPU image' section and repeat/adjust
%%%%%%%%%%% smoothing and threshold parameters until segmentation is
%%%%%%%%%%% satisfactory.


%% generate en face of low DOPU region, and also extract values like mean, thickness, an values
disp('Preparing low DOPU en face...');

% size values
numPoints=size(DOPU_test,1);
numAlines=size(DOPU_test,2);
numBscans=size(DOPU_test,3);

% parameters for the low DOPU region
lowThickness    = zeros(numAlines,numBscans);
lowMean         = zeros(numAlines,numBscans);
totValues       = ones(size(DOPU_test));

for B=1:numBscans % loop across the frames
    for A=1:numAlines % loop across the A-lines
        
        % extract the indices (and round them)
        top = round(DOPU_depth_C(A,B));
        bot = round(DOPU_depth_C_bot(A,B));
        
        % extract the low DOPU region of interest
        lowROI = DOPU_test(top:bot,A,B);
        
        lowThickness(A,B) = length(lowROI); %thickness of the region
        lowMean(A,B) = mean(lowROI); % mean value in that region
        totValues(top:bot,A,B) = lowROI; % low DOPU region values
    end
    fprintf('DOPU frame: %d\n', B);
end

lowValues = totValues(totValues~=1); % save only the low values (space independent)

% display the thickness and mean value maps
figure;
% imagesc(rot90(imadjust(mat2gray(lowThickness,[0,40]),[0,1]),3),[0,1]);colormap(gray);colorbar;axis equal;axis tight;
imagesc(imadjust(mat2gray(lowThickness,[0, 90]),[0,1]),[0,1]);colormap(gray);colorbar;axis equal;axis tight;

title('Thickness');
% ACTUAL max is 50, used 40 for visual clarity (This was for mice!)

figure;
% imagesc(rot90(mat2gray(lowMean,[0,1]),3),[0,1]);colormap(cmap_RPE_r);colorbar;axis equal; axis tight;
imagesc(mat2gray(lowMean,[0,1]),[0,1]);colormap(cmap_RPE_r);colorbar;axis equal; axis tight;

% display with yellow as high dopu and blue/magenta as low
title('DOPU value map');


% save the parameters
save(fullfile(savepath,[filename,'_lowThickness.mat']), 'lowThickness', '-v7.3');
save(fullfile(savepath,[filename,'_lowMean.mat']), 'lowMean', '-v7.3');
save(fullfile(savepath,[filename,'_lowValues.mat']), 'lowValues', '-v7.3');


%% extract OCTA 
% determine the right segmentation area/depths using the 'cehck DOPU image'
% section (the limits are set there). set the top crop startZ to abut 40-100 or less to remove FPN

testOCTA = flipud(OCTA);
startZ=120;
% botlim = 5;
lowlim = 13;
uplim = 20;

fund_all    = zeros(size(cplxOCTA,2),size(cplxOCTA,3));
fund_sup    = zeros(size(cplxOCTA,2),size(cplxOCTA,3));
fund_deep   = zeros(size(cplxOCTA,2),size(cplxOCTA,3));

ind_all     = zeros(size(cplxOCTA,2),size(cplxOCTA,3));
ind_sup     = zeros(size(cplxOCTA,2),size(cplxOCTA,3));
ind_deep    = zeros(size(cplxOCTA,2),size(cplxOCTA,3));

for i = 1:size(cplxOCTA,2)
    for j = 1:size(cplxOCTA,3)
        % extract the indice (and round it)
        top = round(DOPU_depth_C(i,j));
        
%         [fund_all(i,j),ia] = max(cplxOCTA(startZ:top-lowlim,i,j),[],1); 
%         [fund_sup(i,j),is] = max(cplxOCTA(startZ:top-uplim,i,j),[],1);
%         [fund_deep(i,j),id] = max(cplxOCTA(top-uplim:top-lowlim,i,j),[],1);
        
        [fund_all(i,j),ia] = max(testOCTA(startZ:top-lowlim,i,j),[],1); 
        [fund_sup(i,j),is] = max(testOCTA(startZ:top-uplim,i,j),[],1);
        [fund_deep(i,j),id] = max(testOCTA(top-uplim:top-lowlim,i,j),[],1);
        
        % get indices as distances from the top of the RPE
        ind_all(i,j) = top-(ia+startZ-1);
        ind_sup(i,j) = top-(is+startZ-1);
        ind_deep(i,j) = top-(id+top-uplim-1);
    end
end

figure;

OCTA_ALL = rot90(notchfilter(imcenhance(fund_all)),3);
OCTA_ALL_r = imadjust(OCTA_ALL,[0,0.5],[]);
subplot(2,2,1);imagesc(OCTA_ALL_r);colormap(gray);axis equal;
title('OCTA All');

OCTA_deep = rot90(notchfilter(imcenhance(fund_deep)),3); 
OCTA_deep_r = imadjust(OCTA_deep,[0,0.5],[]);
subplot(2,2,2);imagesc(OCTA_deep_r);colormap(gray);axis equal;
title('OCTA deep');

OCTA_sup = rot90(notchfilter(imcenhance(fund_sup)),3); 
OCTA_sup_r = imadjust(OCTA_sup,[0,0.5],[]);
subplot(2,2,3);imagesc(OCTA_sup_r);colormap(gray);axis equal;
title('OCTA sup');

% combine the all (or sup) and dep for higher contrast
fund_comb = max(OCTA_ALL_r,OCTA_deep_r);
OCTA_comb = notchfilter(imcenhance(fund_comb)); 
OCTA_comb_r = imadjust(OCTA_comb,[0,0.5],[]); 
subplot(2,2,4);imagesc(OCTA_comb_r);colormap(gray);axis equal;
title('OCTA comb');


save(fullfile(savepath,[filename,'_OCTA_comb_r.mat']), 'OCTA_comb_r', '-v7.3');
imwrite(OCTA_comb_r,fullfile(savepath,[filename,'_OCTA_comb_r.tif']));

figure;imagesc(rot90(ind_all,3),[0,55]);colormap(cmap_OCTA_r2);axis equal;colorbar;


save(fullfile(savepath,[filename,'_ind_all.mat']), 'ind_all', '-v7.3');


%% combine en face PAF

% combine using the lowThickness as the alpha channel;
figure('position' , [150 150 600 500])
axBG = axes('Position',[0 0 1 1]);

imshow(zeros(numAlines,numBscans)); % prepare a black background

% preparing the mean value map with colourmap of choice
colorMapForSaving = imresize(cmap_RPE_r,[256, 3], 'nearest');
imgM = uint8(255.0 .*mat2gray(lowMean,[0,1]));
imgM_rgb = ind2rgb(imgM,colorMapForSaving);
hold on;
axRPE=axes('Position',[0 0 1 1]);
m=imshow(imgM_rgb); % display the mean value map
% m=imagesc(mat2gray(lowMean,[0,1]),[0,1]); colormap(axRPE,cmap_RPE_r);% display the mean value map

hold off;
set(m, 'AlphaData', imadjust(mat2gray(lowThickness,[0,90]),[0,1])); 
% mask the intensity using the thickness map
% set to 40 because the desired max of 50 is actually too dark to see
% properly.  Should realize that this is indeed 50 max instead
% CHANGE FOR HUMAN DATA TO 90)

print(fullfile(savepath,[filename,'_PAF_DOPU.tif']),'-dtiffn') % save the tiff



% grab next axis and link together
axOCTA = axes('Position',[0 0 1 1]);
linkaxes([axBG,axRPE,axOCTA])
colormap(axRPE , cmap_RPE_r);
colormap(axBG , gray);

colorMapForSaving_r = (imresize(cmap_OCTA_r2,[256, 3], 'nearest'));
colormap(axOCTA , colorMapForSaving_r);

% add vessels from OCTA
hold on;
imgR = uint8(255.0 .*mat2gray(ind_all,[0,55]));
imgR_rgb = ind2rgb(imgR,colorMapForSaving_r);
% axOCTA=axes('Position',[0 0 1 1]);
r=imshow(imgR_rgb);
hold off;
mask = imresize(OCTA_comb_r,[numAlines,numBscans]); % mask intensity with OCTA
set(r, 'AlphaData', mask);




% show colourbars, labels, etc.
set([axBG,axRPE,axOCTA],'Position',[.17 .15 .67 .815]); %from left, from botom, witdh, height
Axis1 = colorbar(axRPE,'Position',[.1 .16 .06 .8]);
Axis2 = colorbar(axBG,'Location','southoutside','Position',[.17 .09 .67 .06]);
Axis3 = colorbar(axOCTA,'Position',[.85 .16 .06 .8]);

Axis1.Label.String= 'DOPU value';

Axis2.Ticks = linspace(0, 1, 11) ; %Create 11 ticks from min to max thickness
labels2 = round(linspace(0,50,11));



Axis2.TickLabels = num2cell(labels2);
Axis2.Label.String= 'Thickness of RPE/choroid in pixels (intensity of background)';

Axis3.Ticks = linspace(0, 1, 12) ; %Create 12 ticks from min to max depth
labels3 = round(linspace(0,55,12));
Axis3.TickLabels = num2cell(labels3);
Axis3.Label.String= 'Vessel distance from RPE in pixels';

title(filename,'Interpreter','none');

print(fullfile(savepath,[filename,'_PAF_full.tif']),'-dtiffn') % save full PAF tiff





%% histogram analysis on the en face low DOPU region
types = {'Thickness','DOPU values'};

for t = 1:2
    hist_type = types{t};

    if strcmp(hist_type,'Thickness')
        value = lowThickness(:);
        binlimits = [0,60];
        number = 30;
        textbox = [0.65 0.7 0.25 0.22];
    else
        value = lowValues(:);
        binlimits = [0,1];
        number = 100;
        textbox = [0.15 0.7 0.25 0.22];
    end

    % title_text = sprintf('Thickness & %s ' ,filename);
    title_text = strcat(hist_type,', volume: ',filename);


    figure;
    h_plot = histogram(value,'Normalization','probability','BinLimits', binlimits, 'NumBins',number); % set bin limits for thickness (around 80?) so the legend isn't wonky


    % Calculate the min, max, mean, median, and standard deviation
    mn=min(value);
    mx=max(value);
    me=nanmean(value);
    md=nanmedian(value);
    stdv=nanstd(value);

    % Create the labels
    minlabel=sprintf('Min - %3.2f', mn);
    maxlabel=sprintf('Max - %3.2f', mx);
    mnlabel=sprintf('Mean - %3.2f', me);
    mdlabel=sprintf('Median - %3.2f', md);
    stdlabel=sprintf('Std Deviation - %3.2f', stdv);
    % Create the textbox
    h=annotation('textbox',textbox); %distace to left, distance to top.  Adjust accordingly
    set(h,'String',{minlabel, maxlabel,mnlabel, mdlabel, stdlabel});
    title(title_text,'Interpreter', 'none');

    saveas(gcf,fullfile(savepath,[filename,'_hist_',hist_type,'.tif']));
    pause(1)
end



