% Melanoma processing from ProcdData data %
close all; clearvars; clc;
loadloc = 'G:\VCSEL_PDOCT\2022.12.14\MIS';    % directory containing the volume(s)
addpath(genpath('C:\Users\coil_\OneDrive\Documents\Github\yusi-miao\Software-MatOctPipeline\PDOCT_400VCSEL'));
addpath(loadloc);
savepath = loadloc;                 % directory to save the files to
fileIdx = 2;

%%% Preset parameter %%%
mcorr       = 0;
matfile     = 0;     % File type of cplxData
numBMscans  = 1;     % number of BM scans
adapt       = 0;     % 1 (adaptive averaging kernel) / 0 (rigid kernel)
kernel      = [3,5]; % averaging kernel dimensions (axial, lateral)
movie       = 1;     % save movie

%%%%% Find filenames of RAW files to be processed %%%%%
cd(loadloc);
files   = (dir('*.unp'));
fnames  = {files.name}';
[~,filename,~] = fileparts(fnames{fileIdx});
disp(['File name: ' filename]);

%%% Load data %%%
load('cmap_dopu_r.mat');    %load the DOPU colormap
load('cmap_RPE_r.mat');     % load the RPE colourmap for en face projection of low DOPU values
% load(fullfile(colormap_dir,'cmap_OCTA_r2.mat'));   % load the depth-encoded OCTA map
if adapt == 1 % need a prior depth profile generated from rigid kernel DOPU first
    load(fullfile(savepath,[filename,'_DOPU_depth_C.mat']));
    load(fullfile(savepath,[filename,'_DOPU_depth_C_bot.mat']));
%     % Need to have the values be flipped upside-down since the volumes are upside-down
    DOPU_depth_C = numPoints-DOPU_depth_C;
    adapt_angles = zeros(numAlines,numBscans);
    adapt_sizes = zeros(numAlines,numBscans);
end

% Load the cplxOCT file(s) %
disp('Loading files...');
try
    % MATLAB data % 
    load(fullfile(loadloc,[filename,'A.mat']));
    load(fullfile(loadloc,[filename,'B.mat']));
catch
    % PyOCT data % 
    load(fullfile(loadloc,[filename,'pdoctA_r.mat']));
    load(fullfile(loadloc,[filename,'pdoctB_r.mat']));
    load(fullfile(loadloc,[filename,'pdoctA_i.mat']));
    load(fullfile(loadloc,[filename,'pdoctB_i.mat']));
    cplxData_A = double(permute(ProcDataA_r+1j*ProcDataA_i,[1,2,3]));
    cplxData_B = double(permute(ProcDataB_r+1j*ProcDataB_i,[1,2,3]));
end

%% Reference frame %%
% Noise estimation from reference frame %
ref_noiseFrame = 990;
nPos = 50;             % axial position of noise estimation
nWid = 10;              % WIDTH of noise estimation

figure(1),subplot(2,1,1),imagesc(imadjust(mat2gray(20*log10(abs(cplxData_A(:,:,ref_noiseFrame)))))),colormap('gray')
subplot(2,1,2),plot(imadjust(mat2gray(mean(20*log10(abs(cplxData_A(:,:,ref_noiseFrame))),2))))

%% Noise-immune DOPU calculation %%
% Calculate Noise %
refOCT_P = cplxData_A(:,:,ref_noiseFrame);
refOCT_S = cplxData_B(:,:,ref_noiseFrame);
refPhaseOff = repmat(angle(sum(refOCT_S.*conj(refOCT_P),1)), [size(refOCT_S,1) 1]); % new paper; get the value from is/os only? check later
refOCT_Cplx = refOCT_P + (refOCT_S.*exp(-1j.*refPhaseOff));

for i = 1:size(refOCT_P,2)
    nOCT_P_real = var(real(refOCT_P(nPos-nWid:nPos+nWid,i))); % orig 1:20
    nOCT_P_imag = var(imag(refOCT_P(nPos-nWid:nPos+nWid,i)));
    nOCT_P_cplx(i) = nOCT_P_real + nOCT_P_imag; % Should this be actualy complex noise?????
    nOCT_S_real = var(real(refOCT_S(nPos-nWid:nPos+nWid,i)));
    nOCT_S_imag = var(imag(refOCT_S(nPos-nWid:nPos+nWid,i)));
    nOCT_S_cplx(i) = nOCT_S_real + nOCT_S_imag;
end
OCT_PN = median(nOCT_P_cplx)+std(nOCT_P_cplx);
OCT_SN = median(nOCT_S_cplx)+std(nOCT_S_cplx);
noise_PS = [OCT_PN, OCT_SN]; 
med_PS=[median(nOCT_P_cplx) median(nOCT_S_cplx)];
std_PS=[std(nOCT_P_cplx) std(nOCT_S_cplx)];
save(fullfile(savepath,[filename,'_noise.mat']), 'noise_PS', 'med_PS', 'std_PS','-v7.3');

% Motion correction %
if mcorr == 1 
    [~, yShift_global, ~] = MotionCorrectionGlobal(cplxData_A(30:end,:,:), filename, savepath);
    for I=1:size(ProcdData_ChP,3)
        cplxData_A(:,  :, I) = circshift(cplxData_A(:, :, I), [yShift_global(I), 0]);
        cplxData_B(:,  :, I) = circshift(cplxData_B(:, :, I), [yShift_global(I), 0]);
    end
end
% figure(2);plot(yShift_global);pause(1);

%% Cropping %%
figure(2);
depthROI = [1 500];
for ii = 1:25:size(cplxData_A,3)
    imagesc(imadjust(mat2gray(20*log10(abs(cplxData_A(:,:,ii)))))),colormap('gray'), hold on
    plot(depthROI(1)*ones(1,size(cplxData_A,2)),'y','LineWidth',3),
    plot(depthROI(2)*ones(1,size(cplxData_A,2)),'y','LineWidth',3),
    pause(0.1)
end

%% DOPU volume processing %%
try
    volume_ChP = cplxData_A(depthROI(1):depthROI(2),:,:);
    volume_ChS = cplxData_B(depthROI(1):depthROI(2),:,:);
catch
    display('No cplxData found')
end 
[numPoints,numAlines,numBscans] = size(volume_ChP);
clearvars cplxData_A cplxData_B

iF=1; % frame counter
clearvars dopu CompCplx
for I = 1:numBMscans:numBscans
    
    K = ((I-1)/numBMscans)+1;
    % Jones Matrix Calculus % 
    for J = 1:numBMscans
        OCT_P  = volume_ChP(:,:,(I+J-1));  
        OCT_S  = volume_ChS(:,:,(I+J-1));

        % Noise-error-corrected Stokes parameters %
        % OL 2014 "Degree of polarization uniformity with high noise
        % immunity using polarization-sensitive oct"
        S0    = OCT_P.*conj(OCT_P) + OCT_S.*conj(OCT_S);    % I
        S1    = OCT_P.*conj(OCT_P) - OCT_S.*conj(OCT_S);    % Q
        S2    = 2.*real(OCT_P.*conj(OCT_S));                % U
        S3    = 2.*imag(OCT_P.*conj(OCT_S));                % V
        S0_NC = S0 - (OCT_PN + OCT_SN);
        S1_NC = S1 - (OCT_PN - OCT_SN);

        % Spatial kernel % 
        if adapt == 1
            S=cat(3,S0_NC,S1_NC,S2,S3);
            [mS, adapt_angles(:,K), adapt_sizes(:,K)]   =  adaptKernelSmoothing_c(S, kernel,DOPU_depth_C(:,K));
            mS0=mS(:,:,1);
            mS1=mS(:,:,2);
            mS2=mS(:,:,3);
            mS3=mS(:,:,4);
        else
            mS0   =  smooth2DFilter(S0_NC, kernel);
            mS1   =  smooth2DFilter(S1_NC, kernel);
            mS2   =  smooth2DFilter(S2, kernel);
            mS3   =  smooth2DFilter(S3, kernel);
        end
        dopu_Numer  = sqrt(mS1.^2 + mS2.^2 + mS3.^2);
        dopu_Denom  = mS0;
        dopu(:,:,J) = dopu_Numer./dopu_Denom; 

        % Bulk-phase correction %
        rPhaseOff = repmat(angle(sum(OCT_S.*conj(OCT_P),1)), [size(OCT_S,1) 1]);
        OCT_Cplx = (OCT_P) + (OCT_S.*exp(-1j.*rPhaseOff)); % phase supposed to be corrected before Stokes Vector??
%             Xconj    = OCT_Cplx.*conj(refOCT_Cplx);
%             bPhaseOff  = repmat(angle(sum(Xconj,1)), [size(Xconj,1) 1]);
%             CompCplx(:,:,J) = OCT_Cplx.*exp(-1j*bPhaseOff);
        CompCplx(:,:,J) = OCT_Cplx;
    end
    
    % BM scan average %
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
    mDOPU(mDOPU<0) = 1; % noise
    DOPU(:,:,K)   = mDOPU;

    if mod(iF,10) == 0
        fprintf('DOPU main process : %d\n', iF);
    end
    iF = iF+1;
end


%% Post-processing DOPU images %%
% frame = 100;
% Idopu = DOPU(:,:,frame);
% mDOPU          = Idopu;
% mDOPU(mDOPU<0)          = 1;
% % mDOPU          = 1./mDOPU;
% % mDOPU(mDOPU<1) = 1;
% % mDOPU   = 1./mDOPU;
% DOPU_filt=mDOPU;
% DOPU_filt(DOPU_filt>0.95) = 1; % threshold the DOPU
% DOPU_filt2=medfilt2(DOPU_filt,[3 5]);
% figure(1),
% subplot(1,3,1),imagesc(mDOPU,[0,1]);colormap(cmap_dopu_r);
% subplot(1,3,2),imagesc(DOPU_filt,[0,1]);colormap(cmap_dopu_r);
% subplot(1,3,3),imagesc(DOPU_filt2,[0,1]);colormap(cmap_dopu_r);

disp('Filtering DOPU...');
DOPU_filt=DOPU;
DOPU_filt(DOPU_filt>0.95) = 1; % threshold the DOPU
DOPU_test=medfilt3(DOPU_filt,[3 5 3]); %[3 5(or 3) 3] (height, wdith, depth) with B-scans, [3 5 13] without B-scans
imagesc((DOPU_test(:,:,100)),[0,1]); colormap('jet')%colormap(cmap_dopu_r);

% Generate depth profiles %
if adapt == 0
	[DOPU_depth_C,DOPU_depth_C_bot]= genDOPU_depth_C(flipud(DOPU_test));
end

% generate composite DOPU+avgOCT logscale B-scans
DOPU_Bscans = genDOPU_combinedBscans(DOPU_test,avgOCT,cmap_dopu_r);

% Save Mat files %
disp('Saving files...');
% avgOCT %
save(fullfile(savepath,[filename,'_avgOCT.mat']), 'avgOCT', '-v7.3');
exportTiff(flipud(avgOCT), fullfile(savepath,[filename,'_avgOCT']))
% DOPU %
if adapt == 1
    save(fullfile(savepath,[filename,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'_adaptSmooth_dopu.mat']), 'DOPU_test', '-v7.3');
    save(fullfile(savepath,[filename,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'_adaptSmooth_dopuBscan.mat']), 'DOPU_Bscans', '-v7.3');
    exportTiff(DOPU_Bscans, fullfile(savepath,[filename,'_adaptive_DOPU_Bscans']))   
else
    save(fullfile(savepath,[filename,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'_dopu.mat']), 'DOPU_test', '-v7.3');
    save(fullfile(savepath,[filename,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'_dopuBscan.mat']), 'DOPU_Bscans', '-v7.3');
    exportTiff(DOPU_Bscans, fullfile(savepath,[filename,'_DOPU_Bscans']))  
    save(fullfile(savepath,[filename,'_DOPU_depth_C.mat']), 'DOPU_depth_C', '-v7.3');
    save(fullfile(savepath,[filename,'_DOPU_depth_C_bot.mat']), 'DOPU_depth_C_bot', '-v7.3');
end
% En face OCT % 
f=(squeeze(mean(avgOCT,1)));
f(f==0)=30;
imgF=10*log10(f); 
imgF=mat2gray(imgF);
imgF_test=imadjust(imgF,[0.25,0.9],[]);
figure(3);imshow(imgF);colormap(gray);axis equal; axis tight;
imwrite(imgF,fullfile(savepath,[filename,'_EnFaceOCT.tif']));

% Save Montage_movie
if movie == 1
    close all
    disp('Saving movies...');
    Montage_movie_func(loadloc, filename, avgOCT, 250);
end
%% Extra: Display DOPU results %%
close all
figure(1);
for ii = 1:10:size(DOPU_Bscans,4)
    
%     imgD = DOPU_test(:,:,ii);
    imgD_rgb_r = DOPU_Bscans(:,:,:,ii);
    [imgD_ind, cmap]=rgb2ind(imgD_rgb_r,32);
%     imagesc(imgD,[0,1]);colormap(cmap_dopu_r);colorbar;
    imagesc(imgD_ind);colormap(cmap);
    
    if adapt==1
        title(['DOPU, adaptive, [',num2str(kernel(1)),',',num2str(kernel(2)),']']);
    else
        title(['DOPU, rigid, [',num2str(kernel(1)),',',num2str(kernel(2)),']']);
    end
    xlabel(ii);
    pause(0.1)
end

frame = 600;
imgD = DOPU_test(:,:,frame);
imgD_rgb_r = DOPU_Bscans(:,:,:,frame);
[imgD_ind, cmap]=rgb2ind(imgD_rgb_r,32);
figure(2),imshow(imgD,[0,1]);colormap(cmap_dopu_r);colorbar;
figure(3),imshow(imgD_ind);colormap(cmap);

