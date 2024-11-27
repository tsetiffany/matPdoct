close all; clear all; clc;

% loadloc = 'D:\DATA\PS-OCT';
% script_dir = 'D:\MJ\Dropbox\ProgramScripts\MatlabScripts\OCTViewer_Project';

% loadloc = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data\MS_DOPU\H001_OS_lM_PS';
% loadloc = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data\OCTA\J001_macula_newdetection';
% loadloc = 'D:\Temporary Group meeting file\A001';
% loadloc = 'E:\DNN_data\3x3mm';
% loadloc = 'E:\CageTyr_ILM01_LFOV';
% loadloc = 'E:\CageRPE65_ILM_LFOV';
% loadloc = 'E:\Cage_WT_ILM_LFOV';

loadloc = 'E:\MFOV WT NonONHRegion';
% currently using 11_20_21 (L), 11_19_34 (flat), and 11_21_16 (R)
% name of collected same frames (L=400,401, F=400,401, R=475,476) = 11_LFR_Cage_WT_RPE_MFOV_Upper

script_dir = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\JonesMatrixCode_Des_Ver6_ProcForPSSAO-OCT';
addpath(script_dir);

%%%% Set parameters %%%%
numBMscans = 2; %4 for all but wide
numFrames = 1000; %1000 for wide, 2000 for others
strtFrame = 1;
kernel = [3,5]; %orig [3,5] for non WFOV, recent [1,8] & [3,10] for wide FOV (1,5 better)

%%% PS mode or not %%%
PS = 1; % 1 (PS mode) / 0 (Non-PS mode)

%% Load the files
fn_num = '11_21_16-Cage_WT_RPE_MFOV_Upper_TiltR'; %%%%%%%%%%%%% Change depending on file name %%%%%%%%%%%%

if PS == 1  
    fn_ChP = [fn_num, '    _mcorr_ChP_Global'];
    load(fullfile(loadloc,fn_ChP));
%     cplxVol_ChP = volume_mcorr_ChP;
    cplxVol_ChP = volume_mcorr_ChP(:,51:end-50,:);

    clear volume_mcorr_ChP
    
    disp('P-channel mcorr loaded');
    
    fn_ChS = [fn_num,'    _mcorr_ChS_Global'];
    load(fullfile(loadloc,fn_ChS));
%     cplxVol_ChS = volume_mcorr_ChS;
    cplxVol_ChS = volume_mcorr_ChS(:,51:end-50,:);
    clear volume_mcorr_ChS

    disp('S-channel mcorr loaded');
else
    fn = [fn_num,'    _mcorr'];
    load(fullfile(loadloc,fn));
    cplxVol  = volume_mcorr;
    clear volume_mcorr
    
    disp('mcorr loaded');
end



%% Process the S & P mcorr files (3D DOPU)
tic;
iF=1; % frame counter
for I = 1:numBMscans:numFrames
    K = ((I-1)/numBMscans)+1;
    if PS == 1 
        refOCT_P = cplxVol_ChP(:,:,I);
        refOCT_S = cplxVol_ChS(:,:,I);
        refPhaseOff = repmat(angle(sum(refOCT_S.*conj(refOCT_P),1)), [size(refOCT_S,1) 1]);
        refOCT_Cplx = refOCT_P + (refOCT_S.*exp(-1j.*refPhaseOff));
        for i = 1:size(refOCT_P,2)
            nOCT_P_real = var(real(refOCT_P(1:20,i))); % orig 1:20
            nOCT_P_imag = var(imag(refOCT_P(1:20,i)));
            nOCT_P_cplx(i) = nOCT_P_real + nOCT_P_imag;
            nOCT_S_real = var(real(refOCT_S(1:20,i)));
            nOCT_S_imag = var(imag(refOCT_S(1:20,i)));
            nOCT_S_cplx(i) = nOCT_S_real + nOCT_S_imag;
        end
        OCT_PN = median(nOCT_P_cplx);
        OCT_SN = median(nOCT_S_cplx);
        for J = 1:numBMscans
            OCT_P  = cplxVol_ChP(:,:,(I+J-1));  
            OCT_S  = cplxVol_ChS(:,:,(I+J-1));
            S0    = OCT_P.*conj(OCT_P) + OCT_S.*conj(OCT_S);
            S1    = OCT_P.*conj(OCT_P) - OCT_S.*conj(OCT_S);
            S2    = 2.*real(OCT_P.*conj(OCT_S));
            S3    = 2.*imag(OCT_P.*conj(OCT_S));
            S0_NC = S0 - (OCT_PN + OCT_SN);
            S1_NC = S1 - (OCT_PN - OCT_SN);
%             figure;subplot(2,2,1);imagesc(S0_NC);colorbar;title('S0_NC');subplot(2,2,2);imagesc(S1_NC);colorbar;title('S1_NC');subplot(2,2,3);imagesc(S2);colorbar;title('S2');subplot(2,2,4);imagesc(S3);colorbar;title('S3');

%             mS0   =  smooth2DFilter(S0_NC, kernel);
%             mS1   =  smooth2DFilter(S1_NC, kernel);
%             mS2   =  smooth2DFilter(S2, kernel);
%             mS3   =  smooth2DFilter(S3, kernel);
            
            mS0   =  adaptKernelSmoothing(S0_NC, kernel,5,0);
            mS1   =  adaptKernelSmoothing(S1_NC, kernel,5,0);
            mS2   =  adaptKernelSmoothing(S2, kernel,5,0);
            mS3   =  adaptKernelSmoothing(S3, kernel,5,0);
%             figure;subplot(2,2,1);imagesc(mS0);colorbar;title('mS0');subplot(2,2,2);imagesc(mS1);colorbar;title('mS1');subplot(2,2,3);imagesc(mS2);colorbar;title('mS2');subplot(2,2,4);imagesc(mS3);colorbar;title('mS3');

            dopu_Numer  = sqrt(mS1.^2 + mS2.^2 + mS3.^2);
            dopu_Denom  = mS0;
            dopu(:,:,J) = dopu_Numer./dopu_Denom; 
            rPhaseOff = repmat(angle(sum(OCT_S.*conj(OCT_P),1)), [size(OCT_S,1) 1]);
            OCT_Cplx = (OCT_P) + (OCT_S.*exp(-1j.*rPhaseOff));
            Xconj    = OCT_Cplx.*conj(refOCT_Cplx);
            bPhaseOff  = repmat(angle(sum(Xconj,1)), [size(Xconj,1) 1]);
            CompCplx(:,:,J) = OCT_Cplx.*exp(-1j*bPhaseOff);
        end
        
        
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
        
        mDOPU          = 1./mDOPU;
        mDOPU(mDOPU<1) = 1;
        DOPU(:,:,K)    = 1./mDOPU;
        
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

%% DOPU serparable 2D space-time kernel
% 2D along BM-scans and time instead. Kernel sizes separate, only within one set of 4-bm scans, no spatial overlap
tic;
kernel = [2,5]; % no smoothing along depth this time
iF=1; % frame counter
for I = 1:numBMscans:numFrames
    K = ((I-1)/numBMscans)+1;
    if PS == 1 
        refOCT_P = cplxVol_ChP(:,:,I);
        refOCT_S = cplxVol_ChS(:,:,I);
        refPhaseOff = repmat(angle(sum(refOCT_S.*conj(refOCT_P),1)), [size(refOCT_S,1) 1]);
        refOCT_Cplx = refOCT_P + (refOCT_S.*exp(-1j.*refPhaseOff));
        for i = 1:size(refOCT_P,2)
            nOCT_P_real = var(real(refOCT_P(1:20,i))); % orig 1:20
            nOCT_P_imag = var(imag(refOCT_P(1:20,i)));
            nOCT_P_cplx(i) = nOCT_P_real + nOCT_P_imag;
            nOCT_S_real = var(real(refOCT_S(1:20,i)));
            nOCT_S_imag = var(imag(refOCT_S(1:20,i)));
            nOCT_S_cplx(i) = nOCT_S_real + nOCT_S_imag;
        end
        OCT_PN = median(nOCT_P_cplx);
        OCT_SN = median(nOCT_S_cplx);
        
        %smooth across width of each BM scan
        for J = 1:numBMscans
            OCT_P  = cplxVol_ChP(:,:,(I+J-1));  
            OCT_S  = cplxVol_ChS(:,:,(I+J-1));
            S0    = OCT_P.*conj(OCT_P) + OCT_S.*conj(OCT_S);
            S1    = OCT_P.*conj(OCT_P) - OCT_S.*conj(OCT_S);
            S2    = 2.*real(OCT_P.*conj(OCT_S));
            S3    = 2.*imag(OCT_P.*conj(OCT_S));
            S0_NC = S0 - (OCT_PN + OCT_SN);
            S1_NC = S1 - (OCT_PN - OCT_SN);
%             figure;subplot(2,2,1);imagesc(S0_NC);colorbar;title('S0_NC');subplot(2,2,2);imagesc(S1_NC);colorbar;title('S1_NC');subplot(2,2,3);imagesc(S2);colorbar;title('S2');subplot(2,2,4);imagesc(S3);colorbar;title('S3');

%             mS0   =  smooth2DFilter(S0_NC, kernel);
%             mS1   =  smooth2DFilter(S1_NC, kernel);
%             mS2   =  smooth2DFilter(S2, kernel);
%             mS3   =  smooth2DFilter(S3, kernel);
            
            mS0   =  adaptKernelSmoothing(S0_NC, kernel,5,0);
            mS1   =  adaptKernelSmoothing(S1_NC, kernel,5,0);
            mS2   =  adaptKernelSmoothing(S2, kernel,5,0);
            mS3   =  adaptKernelSmoothing(S3, kernel,5,0);
%             figure;subplot(2,2,1);imagesc(mS0);colorbar;title('mS0');subplot(2,2,2);imagesc(mS1);colorbar;title('mS1');subplot(2,2,3);imagesc(mS2);colorbar;title('mS2');subplot(2,2,4);imagesc(mS3);colorbar;title('mS3');

            dopu_Numer  = sqrt(mS1.^2 + mS2.^2 + mS3.^2);
            dopu_Denom  = mS0;
            dopu(:,:,J) = dopu_Numer./dopu_Denom; 
            rPhaseOff = repmat(angle(sum(OCT_S.*conj(OCT_P),1)), [size(OCT_S,1) 1]);
            OCT_Cplx = (OCT_P) + (OCT_S.*exp(-1j.*rPhaseOff));
            Xconj    = OCT_Cplx.*conj(refOCT_Cplx);
            bPhaseOff  = repmat(angle(sum(Xconj,1)), [size(Xconj,1) 1]);
            CompCplx(:,:,J) = OCT_Cplx.*exp(-1j*bPhaseOff);
        end
        
        
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
        
        mDOPU          = 1./mDOPU;
        mDOPU(mDOPU<1) = 1;
        DOPU(:,:,K)    = 1./mDOPU;
        
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


%% DOPU 2D kernel (averaging first followed by DOPU calculation)
kernel = [6,6];
tic;
iF=1; % frame counter
for I = 1:numBMscans:numFrames
    K = ((I-1)/numBMscans)+1;

    % compensate the phase between the two channels
    refOCT_P = cplxVol_ChP(:,:,I);
    refOCT_S = cplxVol_ChS(:,:,I);
%     refPhaseOff = repmat(angle(sum(refOCT_S.*conj(refOCT_P),1)), [size(refOCT_S,1) 1]);
%     refOCT_Cplx = refOCT_P + (refOCT_S.*exp(-1j.*refPhaseOff));

    % noise calculation (from non-averaged OCTs
    for i = 1:size(refOCT_P,2)
        nOCT_P_real = var(real(refOCT_P(1:20,i))); % orig 1:20
        nOCT_P_imag = var(imag(refOCT_P(1:20,i)));
        nOCT_P_cplx(i) = nOCT_P_real + nOCT_P_imag;
        nOCT_S_real = var(real(refOCT_S(1:20,i)));
        nOCT_S_imag = var(imag(refOCT_S(1:20,i)));
        nOCT_S_cplx(i) = nOCT_S_real + nOCT_S_imag;
    end
    OCT_PN = median(nOCT_P_cplx);
    OCT_SN = median(nOCT_S_cplx);
%     
    
    %------------- Method 1 (from old processing code)-----------------
%     % grab the chunk of the volume with all BM scans
%     OCT_P_all = cplxVol_ChP(:,:,I:I+numBMscans-1);
%     OCT_S_all = cplxVol_ChS(:,:,I:I+numBMscans-1);
%     
%     % calculate the BMscan phase offset for each channel and compensate
%     bPhaseOffP = angle(OCT_P_all)-angle(refOCT_P);
%     bPhaseOffS = angle(OCT_S_all)-angle(refOCT_S); % or can do per frame with the complex conjugate of the reference and angle that
%     OCT_P = mean(abs(OCT_P_all).*exp(1i.*(angle(OCT_P_all)-bPhaseOffP)), 3);   % entries of gPhase averaged matrix
%     OCT_S = mean(abs(OCT_S_all).*exp(1i.*(angle(OCT_S_all)-bPhaseOffS)), 3);

    % ------------ Method 2 (from current processing code) --------------
    % loop across the BM scans
    for J = 1:numBMscans
        % grab the frame
        OCT_Cplx_P    = cplxVol_ChP(:,:,(I+J-1));
        % find the angle relative to the reference frame (first BM scan)
        Xconj_P    = OCT_Cplx_P.*conj(refOCT_P);
        bPhaseOff_P   = repmat(angle(sum(Xconj_P,1)), [size(Xconj_P,1) 1]);
        % compensate the offset
        OCT_P_all(:,:,J) = OCT_Cplx_P.*exp(-1j*bPhaseOff_P);
        
        % repeat for other channel
        OCT_Cplx_S    = cplxVol_ChS(:,:,(I+J-1));
        Xconj_S    = OCT_Cplx_S.*conj(refOCT_S);
        bPhaseOff_S   = repmat(angle(sum(Xconj_S,1)), [size(Xconj_S,1) 1]);
        OCT_S_all(:,:,J) = OCT_Cplx_S.*exp(-1j*bPhaseOff_S);
    end
    % sum (MEAN!) across the channels
    OCT_P = mean(OCT_P_all, 3);   
    OCT_S = mean(OCT_S_all, 3);


%     % noise calculation (from averaged OCTs)
%     for i = 1:size(OCT_P,2)
%         nOCT_P_real = var(real(OCT_P(1:20,i))); % orig 1:20
%         nOCT_P_imag = var(imag(OCT_P(1:20,i)));
%         nOCT_P_cplx(i) = nOCT_P_real + nOCT_P_imag;
%         nOCT_S_real = var(real(OCT_S(1:20,i)));
%         nOCT_S_imag = var(imag(OCT_S(1:20,i)));
%         nOCT_S_cplx(i) = nOCT_S_real + nOCT_S_imag;
%     end
%     OCT_PN = median(nOCT_P_cplx);
%     OCT_SN = median(nOCT_S_cplx);
    
    
%     
%     OCT_P_orig  = cplxVol_ChP(:,:,I);  
%     OCT_S_orig  = cplxVol_ChS(:,:,I);
    
            
    % calcualte the stokes vector elements
    S0    = OCT_P.*conj(OCT_P) + OCT_S.*conj(OCT_S);
    S1    = OCT_P.*conj(OCT_P) - OCT_S.*conj(OCT_S);
    S2    = 2.*real(OCT_P.*conj(OCT_S));
    S3    = 2.*imag(OCT_P.*conj(OCT_S));
    S0_NC = S0 - (OCT_PN + OCT_SN);
    S1_NC = S1 - (OCT_PN - OCT_SN);
%             figure;subplot(2,2,1);imagesc(S0_NC);colorbar;title('S0_NC');subplot(2,2,2);imagesc(S1_NC);colorbar;title('S1_NC');subplot(2,2,3);imagesc(S2);colorbar;title('S2');subplot(2,2,4);imagesc(S3);colorbar;title('S3');

%             mS0   =  smooth2DFilter(S0_NC, kernel);
%             mS1   =  smooth2DFilter(S1_NC, kernel);
%             mS2   =  smooth2DFilter(S2, kernel);
%             mS3   =  smooth2DFilter(S3, kernel);

    mS0   =  adaptKernelSmoothing(S0_NC, kernel,5,0);
    mS1   =  adaptKernelSmoothing(S1_NC, kernel,5,0);
    mS2   =  adaptKernelSmoothing(S2, kernel,5,0);
    mS3   =  adaptKernelSmoothing(S3, kernel,5,0);
%             figure;subplot(2,2,1);imagesc(mS0);colorbar;title('mS0');subplot(2,2,2);imagesc(mS1);colorbar;title('mS1');subplot(2,2,3);imagesc(mS2);colorbar;title('mS2');subplot(2,2,4);imagesc(mS3);colorbar;title('mS3');

    dopu_Numer  = sqrt(mS1.^2 + mS2.^2 + mS3.^2);
    dopu_Denom  = mS0;
    mDOPU = dopu_Numer./dopu_Denom; 

    mDOPU          = 1./mDOPU;
    mDOPU(mDOPU<1) = 1;
    DOPU(:,:,K)    = 1./mDOPU;
        
    
    fprintf('Frame processed : %d\n', iF);
    iF = iF+1;
end
toc


%% Save the images
disp('Saving files...');
if PS == 1
%     save(fullfile(loadloc,[fn_num,'_DOPU_',num2str(kernel(1)),'_',num2str(kernel(2)),'.mat']), 'DOPU', '-v7.3');
%     save(fullfile(loadloc,[fn_num,'_DOPU_',num2str(kernel(1)),'_',num2str(kernel(2)),'_adaptSmooth.mat']), 'DOPU', '-v7.3');
%     save(fullfile(loadloc,[fn_num,'_DOPU_',num2str(kernel(1)),'_',num2str(kernel(2)),'_adaptSmooth_2D_r1.mat']), 'DOPU', '-v7.3');
    save(fullfile(loadloc,[fn_num,'_DOPU_',num2str(kernel(1)),'_',num2str(kernel(2)),'_adaptSmooth_2D_TS.mat']), 'DOPU', '-v7.3');
end

% save(fullfile(loadloc,[fn_num,'_avgOCT.mat']), 'avgOCT', '-v7.3');
% save(fullfile(loadloc,[fn_num,'_OCTA.mat']), 'OCTA', '-v7.3');
% save(fullfile(loadloc,[fn_num,'_OCTA_Phs.mat']), 'OCTA_Phs', '-v7.3');
disp('Files saved');


