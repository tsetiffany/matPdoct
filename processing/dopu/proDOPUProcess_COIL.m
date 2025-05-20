close all;
clear all;
clc;
loadloc = 'G:\EDUADO\RP006_2022.05.17\OS';
savepath = loadloc;
fileStart = 1;
depthROI = [1, 550];

proDOPUProcess_COIL_func(loadloc,savepath,fileStart, depthROI);


function proDOPUProcess_COIL_func(loadloc,savepath,fileStart, depthROI)
    %% set up loading directories
    disp('Setting up load directiories...');
    % loadloc = '/project/6007991/borg/STUDENTS/ymiao/Melanoma/MJJ_Logitudianl/2022.02.23_MJ';
    addpath(loadloc);
    disp(loadloc);
    % savepath='/project/6007991/borg/STUDENTS/destinyh/DATA/Melanoma/MJJ_Logitudianal/2022.02.23_MJ';
    disp(savepath);
    if ~exist(savepath, 'dir')
       mkdir(savepath)
    end
    disp('Adding script directory...');
    script_dir = 'C:\Users\coil_\Desktop\Github\Project-PD_OCT\Melanoma\code';
    addpath(script_dir);                % directory containing this script

    %%%%% Find filenames of RAW files to be processed %%%%%
    disp('Finding filenames of RAW files to be processed...');
    cd(loadloc);
    files   = (dir('*.unp'));
    fnames  = {files.name}';

    fprintf('Starting file loop at file start index %d for number of files %d...\n',fileStart,length(fnames));
    for fileIdx = fileStart:length(fnames)
        %% Load the file(s) %%%
        fn = fnames{fileIdx};
        filename= fn(1:end-5);
        fprintf('--------------------File number %d, name %s --------------------\n',fileIdx,filename);
        fn_num=filename;

        disp('Loading files...');
        % Check if it was processed with MATLAB or Python
        if isfile(fullfile(loadloc,[filename,'-cplxOCT_A.mat']))
            % MATLAB: just load directly
            load(fullfile(loadloc,[filename,'-cplxOCT_A.mat']));
            load(fullfile(loadloc,[filename,'-cplxOCT_B.mat']));
        else
            % Python: recombine complex parts
            load(fullfile(loadloc,[filename,'-pdoctA_r.mat']));
            load(fullfile(loadloc,[filename,'-pdoctA_i.mat']));
            cplxData_A = complex(ProcdDataA_r, ProcdDataA_i);
            clear ProcdDataA_r ProcdDataA_i
            
            load(fullfile(loadloc,[filename,'-pdoctB_r.mat']));
            load(fullfile(loadloc,[filename,'-pdoctB_i.mat']));
            cplxData_B = complex(ProcdDataB_r, ProcdDataB_i);
            clear ProcdData_r ProcdDataB_i
        end
        
        % reset sizes
        [numPoints,numAlines,numBscans] = size(cplxData_A);
        
%         ii=509;
%         imgP = mat2gray(squeeze(20.*log10(abs(cplxData_A(:,:,ii)))));
%         figPreMcorr=figure('WindowState', 'maximized'); 
%         imagesc(imadjust(imgP)); colormap(gray);
%         title([filename, ' preMcorr'],'Interpreter','none');
%         xlabel(num2str(ii));
%         saveas(figPreMcorr,fullfile(savepath,[fn_num,'_preMcorr.tif']));

        %% reference frame noise processing
        disp('Reference frame noise processing...');
        ref_Frame=990;
        refOCT_P = cplxData_A(:,:,ref_Frame);
        refOCT_S = cplxData_B(:,:,ref_Frame);
        refPhaseOff = repmat(angle(sum(refOCT_S.*conj(refOCT_P),1)), [size(refOCT_S,1) 1]); % new paper; get the value from is/os only? check later
        refOCT_Cplx = refOCT_P + (refOCT_S.*exp(-1j.*refPhaseOff));
        
        Nroi = 50;
        for i = 1:size(refOCT_P,2)
            nOCT_P_real = var(real(refOCT_P(Nroi-10:Nroi+10,i))); % orig 1:20
            nOCT_P_imag = var(imag(refOCT_P(Nroi-10:Nroi+10,i)));
            nOCT_P_cplx(i) = nOCT_P_real + (nOCT_P_imag); % Should this be actualy complex noise?????
            nOCT_S_real = var(real(refOCT_S(Nroi-10:Nroi+10,i)));
            nOCT_S_imag = var(imag(refOCT_S(Nroi-10:Nroi+10,i)));
            nOCT_S_cplx(i) = nOCT_S_real + (nOCT_S_imag);
        end

        OCT_PN = median(nOCT_P_cplx)+std(nOCT_P_cplx);
        OCT_SN = median(nOCT_S_cplx)+std(nOCT_S_cplx);

        noise_PS=[OCT_PN, OCT_SN];
        med_PS=[median(nOCT_P_cplx) median(nOCT_S_cplx)];
        std_PS=[std(nOCT_P_cplx) std(nOCT_S_cplx)];
        fprintf('Noise levels, OCT_PN = %.2e, OCT_SN = %.2e \n',OCT_PN, OCT_SN);
        save(fullfile(savepath,[fn_num,'_noise.mat']), 'noise_PS', 'med_PS','std_PS','-v7.3');

        %% Motion correction
        disp('Motion correction...');
        % Cropping
        ProcdData_ChP=cplxData_A(31:end,:,:);
        ProcdData_ChS=cplxData_B(31:end,:,:);
        clear cplxData_A cplxData_B

        % reset sizes
        [numPoints,numAlines,numBscans] = size(ProcdData_ChP);
        
%         addpath('/project/6007991/borg/STUDENTS/destinyh/MATLAB/volume_mcorr_code');
        addpath('C:\Users\coil_\Desktop\Github\Project-PD_OCT\Melanoma\code\volume_mcorr_code');
        [~, yShift_global, ~] = MotionCorrectionGlobal(ProcdData_ChP, filename, savepath);

        volume_mcorr_ChP = zeros(size(ProcdData_ChP));
        volume_mcorr_ChS = zeros(size(ProcdData_ChS));

        for I=1:numBscans
            volume_mcorr_ChS(:,  :, I) = circshift(ProcdData_ChS(:, :, I), [yShift_global(I), 0]);
            volume_mcorr_ChP(:,  :, I) = circshift(ProcdData_ChP(:, :, I), [yShift_global(I), 0]);
        end

    %     % mcorr shifts
    %     figShift=figure(2); 
    %     plot(yShift_global);pause(1);
    %     saveas(figShift,fullfile(savepath,[fn_num,'_mcorr_shifts.tif']));

        % central fast scan
        ii=500;
        figFast=figure(3); 
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
        saveas(figFast,fullfile(savepath,[fn_num,'-fig_mcorr_fast.jpg']));


        % Central slow scan
        figSlow=figure(4); 
        subplot(2,2,1)
        imgP = mat2gray(squeeze(20.*log10(abs(ProcdData_ChP(:,ii,:)))));
        imagesc(imadjust(imgP)); colormap(gray);title('ChP');

        subplot(2,2,2)
        imgPM = mat2gray(squeeze(20.*log10(abs(volume_mcorr_ChP(:,ii,:)))));
        imagesc(imadjust(imgPM)); colormap(gray);title('ChP - mcorr');

        subplot(2,2,3)
        imgS = mat2gray(squeeze(20.*log10(abs(ProcdData_ChS(:,ii,:)))));
        imagesc(imadjust(imgS)); colormap(gray);title('ChS');

        subplot(2,2,4)
        imgSM = mat2gray(squeeze(20.*log10(abs(volume_mcorr_ChS(:,ii,:)))));
        imagesc(imadjust(imgSM)); colormap(gray);title('ChS - mcorr');
        saveas(figSlow,fullfile(savepath,[fn_num,'-fig_mcorr_slow.jpg']));
        
        % saving mcorr files
%         disp('Saving mcorrs...');
%         save(fullfile(savepath,[filename,'    _mcorr_ChP_Global.mat']), 'volume_mcorr_ChP', '-v7.3');
%         save(fullfile(savepath,[filename,'    _mcorr_ChS_Global.mat']), 'volume_mcorr_ChS', '-v7.3');
        save(fullfile(savepath,[filename,'    _mcorr_yShift_global.mat']), 'yShift_global', '-v7.3');
        disp('mcorrs saved');

        %% More cropping if necessary
%         %ctop = 1;
%         %cbot = 400;
%         cbot=ctop+399;
        figCrop=figure(5); 
        ii=500;
        imagesc(imadjust(mat2gray(20*log10(abs(volume_mcorr_ChP(:,:,ii)))))),colormap('gray'), hold on
        plot(depthROI(1)*ones(1,size(volume_mcorr_ChP,2)),'y','LineWidth',3),
        plot(depthROI(2)*ones(1,size(volume_mcorr_ChP,2)),'y','LineWidth',3),
        hold off;
        title('Cropping boundaries');
        saveas(figCrop,fullfile(savepath,[fn_num,'-fig_crop_bounds.jpg']));
        
        volume_mcorr_ChP = volume_mcorr_ChP(depthROI(1):depthROI(2),:,:);
        volume_mcorr_ChS = volume_mcorr_ChS(depthROI(1):depthROI(2),:,:);

        

        %% DOPU process
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load colormap %
        colormap_dir = 'C:\Users\coil_\Desktop\Github\Project-PD_OCT\Melanoma\code\JonesMatrixCode\data';
        load(fullfile(colormap_dir,'cmap_dopu_r.mat'));    %load the DOPU colormap
        load(fullfile(colormap_dir,'cmap_RPE_r.mat'));     % load the RPE colourmap for en face projection of low DOPU values
        % load(fullfile(colormap_dir,'cmap_OCTA_r2.mat'));   % load the depth-encoded OCTA map

        %%%% Set parameters %%%%
        numBMscans = 1;     % number of BM scans
        kernel = [3,5];     % averaging kernel dimensions (axial, lateral)
        
        %%%% Size setup %%%%
        [numPoints,numAlines, numFrames] = size(volume_mcorr_ChS);
        numBscans = round(numFrames./numBMscans);
        
        DOPU = zeros(numPoints,numAlines,numBscans);
        avgOCT = zeros(numPoints,numAlines,numBscans);
        OCTA = zeros(numPoints,numAlines,numBscans);

        % % volumes are expected to be upside-down. Flip if they started rightside-up
        % volume_mcorr_ChP = flipud(volume_mcorr_ChP);
        % volume_mcorr_ChS = flipud(volume_mcorr_ChS);

        %% processing the files
        adapt = 0; % do rigid first, then repeat with adaptive
        while adapt <=1
            tic;
            iF=1; % frame counter
            for I = 1:numBMscans:numFrames
                K = ((I-1)/numBMscans)+1;
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
                        % new adaptive kernel with overall depth map
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
                    rPhaseOff = repmat(angle(sum(OCT_S.*conj(OCT_P),1)), [size(OCT_S,1) 1]);
                    OCT_Cplx = (OCT_P) + (OCT_S.*exp(-1j.*rPhaseOff));
%                     Xconj    = OCT_Cplx.*conj(refOCT_Cplx);
%                     bPhaseOff  = repmat(angle(sum(Xconj,1)), [size(Xconj,1) 1]);
%                     CompCplx(:,:,J) = OCT_Cplx.*exp(-1j*bPhaseOff);
                    CompCplx(:,:,J) = OCT_Cplx;
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


                fprintf('Frame processed : %d\n', iF);
                iF = iF+1;
            end
            toc
            clear dopu CompCplx
            
            % Post-processing DOPU images %
            disp('Filtering DOPU...');
            DOPU_filt=DOPU;
            DOPU_filt(DOPU_filt>0.95) = 1; % threshold the DOPU
            DOPU_test=medfilt3(DOPU_filt,[3 5 3]); %[3 5(or 3) 3] (height, wdith, depth) with B-scans, [3 5 13] without B-scans
            disp('DOPU filtered.');
            % imagesc((DOPU_test(:,:,5)),[0,1]);colormap(cmap_dopu_r);
                
            if adapt == 0
                % Save the files for rigid and avg
                disp('Saving DOPU, avgOCT files...');
                % generate composite DOPU+avgOCT logscale B-scans
                DOPU_Bscans = genDOPU_combinedBscans(DOPU_test,avgOCT,cmap_dopu_r);
                save(fullfile(savepath,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'_dopu.mat']), 'DOPU_test', '-v7.3');
                save(fullfile(savepath,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'_dopuBscan.mat']), 'DOPU_Bscans', '-v7.3');
                exportTiff(DOPU_Bscans, fullfile(savepath,[fn_num,'_DOPU_Bscans']))
    
                save(fullfile(savepath,[fn_num,'_avgOCT.mat']), 'avgOCT', '-v7.3');
                exportTiff(flipud(avgOCT), fullfile(savepath,[fn_num,'_avgOCT']))
                save(fullfile(savepath,[fn_num,'_OCTA.mat']), 'OCTA', '-v7.3');
                disp('Files saved');
            
                

                % generate depth profiles
                [DOPU_depth_C,DOPU_depth_C_bot]= genDOPU_depth_C(flipud(DOPU_test));
                save(fullfile(savepath,[filename,'_DOPU_depth_C.mat']), 'DOPU_depth_C', '-v7.3');
                save(fullfile(savepath,[filename,'_DOPU_depth_C_bot.mat']), 'DOPU_depth_C_bot', '-v7.3');
                
                % save an example of the cut lines
                ii=500;
                fig_DOPUdepth=figure;
                imgA = squeeze(flipud(DOPU_test(:,:,ii)));
                imagesc(imgA,[0,1]);colormap(cmap_dopu_r);colorbar;%colormap(gray, cmap_dopu);
                title('DOPU depth estimations');
                xlabel(num2str(ii));
                hold on;
                plot(DOPU_depth_C(:,ii),'black','LineWidth', 2);
                plot(DOPU_depth_C_bot(:,ii),'black','LineWidth', 2);
                hold off;
                saveas(fig_DOPUdepth,fullfile(savepath,[fn_num,'-fig_DOPU_bounds.jpg']));
                
                % Need to have the values be flipped upside-down since the volumes are upside-down
                DOPU_depth_C = numPoints-DOPU_depth_C;
                adapt_angles = zeros(numAlines,numBscans);
                adapt_sizes = zeros(numAlines,numBscans);
                
            elseif adapt == 1
                disp('Saving adaptive DOPU...')
                % generate composite DOPU+avgOCT logscale B-scans
                DOPU_Bscans = genDOPU_combinedBscans(DOPU_test,avgOCT,cmap_dopu_r);
                save(fullfile(savepath,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'_adaptSmooth_dopu.mat']), 'DOPU_test', '-v7.3');
                save(fullfile(savepath,[fn_num,'_DOPU_Bscans_',num2str(kernel(1)),'_',num2str(kernel(2)),'_adaptSmooth_dopuBscan.mat']), 'DOPU_Bscans', '-v7.3');
                exportTiff(DOPU_Bscans, fullfile(savepath,[fn_num,'_adaptive_DOPU_Bscans']))
                disp('Files saved.')
            end
            
            
            
            adapt = adapt +1;
        end

        %%%%%%%%%%%%%  Display the images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save En face OCT image% 
        f=(squeeze(mean(avgOCT,1)));
        f(f==0)=30;
        imgF=10*log10(f); 
        imgF=mat2gray(imgF);
        imgF_test=imadjust(imgF,[0.25,0.9],[]);
        figure(3);imshow(imgF);colormap(gray);axis equal; axis tight;
        imwrite(imgF,fullfile(savepath,[fn_num,'_EnFaceOCT.jpg']));

        % Save Montage_movie
        Montage_movie_func(loadloc, filename, avgOCT);

        fprintf('----------------OCT volume %d %s process complete-----------------\n\r', fileIdx,filename);

    end


end
