function [cplxData_A, cplxData_B] = process_oct_volume(fn,LUT,parameters,dispMaxOrder,dispCoeffs)
    %% Volume process %%
    tic
    clearvars ProcdData

    depthROI        = [1, 1000];
    bitDepth        = 12;
    byteSize        = bitDepth/8;

    %%% Acquisition Parameters %%%
    numPoints     = parameters(1)*2;
    numAscans     = parameters(2)/2;
    numBscans     = parameters(3);
    numCscans     = parameters(4);
    numMscans     = parameters(5);

    fid = fopen(fn);
    ProcdData = zeros(round(numPoints/2),numAscans*2,numBscans);
    counter = 1;
    for FrameNum = 1:numBscans
%         fseek(fid,2*numPoints*numAscans*(FrameNum-1),-1);
%         rawData_interlace = fread(fid,[numPoints,numAscans], 'ubit12');
        if bitDepth == 16
            fid = fopen(fn);
            fseek(fid,byteSize*numPoints*numAscans*(FrameNum-1),-1);
            rawData_interlace = fread(fid,[numPoints,numAscans], 'uint16');
            fclose(fid);
        elseif bitDepth == 12
            rawData_interlace = unpack_u12u16(fn,numPoints,numAscans,FrameNum);
        else
            error('%d-bit data not supported.', bitDepth)
        end
        rawData = hilbert(cat(2, rawData_interlace(1:2:end,:),rawData_interlace(2:2:end,:)));
        
        % Rescale process %
        rawData_Rescaled = reSampling_LUT(rawData,LUT);
%         fftData_Rescaled = fft(rawData_Rescaled);
        
        % FPN remove %
        rawData_FPNSub_A = rawData_Rescaled(:,1:end/2)...
            - (repmat(median(real(rawData_Rescaled(:,1:end/2)),2), [1,size(rawData_Rescaled(:,1:end/2),2)])...
            +1j.*repmat(median(imag(rawData_Rescaled(:,1:end/2)),2), [1,size(rawData_Rescaled(:,1:end/2),2)]));
        rawData_FPNSub_B = rawData_Rescaled(:,end/2+1:end)...
            - (repmat(median(real(rawData_Rescaled(:,end/2+1:end)),2), [1,size(rawData_Rescaled(:,end/2+1:end),2)])...
            +1j.*repmat(median(imag(rawData_Rescaled(:,end/2+1:end)),2), [1,size(rawData_Rescaled(:,end/2+1:end),2)]));
        rawData_FPNSub = cat(2,rawData_FPNSub_A,rawData_FPNSub_B);
        
        % Windowing process %
        rawData_HamWin = rawData_FPNSub...
            .*repmat(hann(size(rawData_FPNSub,1)),[1 size(rawData_FPNSub,2)]);
    
        % Dispersion estimation & compensation %
        rawData_DisComp = compDisPhase(rawData_HamWin,dispMaxOrder,dispCoeffs);
        fftData_DispComp = fft(rawData_DisComp);
    
        ProcdData(:,:,FrameNum) = fftData_DispComp(1:round(numPoints/2),:);
%         fprintf('OCT volume process: %d\n', FrameNum);
        if mod(counter,10) == 0
            fprintf('- OCT volume process: %d /%d\n ', FrameNum, numBscans);
        end
        counter = counter+1;
    end
%     toc
    % save(fullfile(cd,fname_save), 'ProcdData', '-v7.3');
    cplxData_A = reorder_bscan(ProcdData(depthROI(1):depthROI(2),1:end/2,:),numMscans);
    cplxData_B = reorder_bscan(ProcdData(depthROI(1):depthROI(2),end/2+1:end,:),numMscans);
%     cplxData_A = (ProcdData(depthROI(1):depthROI(2),1:end/2,:));
%     cplxData_B = (ProcdData(depthROI(1):depthROI(2),end/2+1:end,:));
    fclose(fid);

%     [~,fname_save,~] = fileparts(fn);
%     exportTiff(flip(20*log10(abs(cplxData_A)+abs(cplxData_B))),fullfile(savepath,[fname_save,'_avgOCT']))
    % exportTiff(flip(20*log10(abs(cplxData_B))),[fname_save,'OCT_B'])
%     save(fullfile(savepath,[fname_save,'A']), 'cplxData_A', '-v7.3');
%     save(fullfile(savepath,[fname_save,'B']), 'cplxData_B', '-v7.3');

end
