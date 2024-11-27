function [avgOCT,DOPU, OCTA] = process_dopu_volume(cplxData_A,cplxData_B, OCT_PN, OCT_SN,numBMscans)
%%% Preset parameter %%%
depthROI    = [1 1000];
adapt       = 0;          % 1 (adaptive averaging kernel) / 0 (rigid kernel)
kernel      = [3,5];     % averaging kernel dimensions (axial, lateral)
filtDOPU    = 1;

    % Crop volume %
    try
        volume_mcorr_ChP = cplxData_A(depthROI(1):depthROI(2),:,:);
        volume_mcorr_ChS = cplxData_B(depthROI(1):depthROI(2),:,:);
    catch
        volume_mcorr_ChP = cplxData_A;
        volume_mcorr_ChS = cplxData_B;
        disp('depthROI out of range: no crop applied')
    end
    clearvars cplxData_A cplxData_B
    % reset sizes
    [numPoints,numAlines,numBscans] = size(volume_mcorr_ChP);

    %%% MAIN DOPU %%%
    %     disp('DOPU Main Process...')
    iF=1; % frame counter
    clearvars dopu CompCplx
    for I = 1:numBMscans:numBscans
        
        K = ((I-1)/numBMscans)+1;
        % Jones Matrix Calculus % 
        for J = 1:numBMscans
	    if I+J-1>numBscans
		break
	    end
            OCT_P  = volume_mcorr_ChP(:,:,(I+J-1));  
            OCT_S  = volume_mcorr_ChS(:,:,(I+J-1));
    
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

    %%% Post-processing DOPU  %%%
    if filtDOPU == 1
        disp('Filtering DOPU...');
        DOPU_filt=DOPU;
        DOPU_filt(DOPU_filt>0.95) = 1; % threshold the DOPU
        DOPU=medfilt3(DOPU_filt,[3 5 3]); %[3 5(or 3) 3] (height, wdith, depth) with B-scans, [3 5 13] without B-scans
    end
    
%     % generate composite DOPU+avgOCT logscale B-scans
%     DOPU_Bscans = genDOPU_combinedBscans(DOPU,avgOCT,cmap_dopu_r);

    disp('DOPU Processing complete');
    close all
end
