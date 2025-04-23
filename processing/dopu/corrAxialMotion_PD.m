function corrAxialMotion_PD(fn,volume_ChP, volume_ChS,folder)
    %% Compute maximum cross-correlation (axial only) 
    motionA = maxxcorrAx((abs(volume_ChP).^2));
    xaxis = 1:1:size(motionA,2);

    %% Set smoothing parameter
    p = polyfit(xaxis,motionA,2);
    f = polyval(p,xaxis);

    %% Compute motion correction parameters and do motion correction       
    disp_ind = motionA - f;     
    m = size(volume_ChP,1);
    n = size(volume_ChP,3);
    topZero = max(disp_ind);
    botZero = abs(min(disp_ind));
    for k=1:n
        top = round(topZero-disp_ind(k));
        volume_mcorr_ChP(top+1:top+m,:,k) = volume_ChP(:,:,k);
        volume_mcorr_ChS(top+1:top+m,:,k) = volume_ChS(:,:,k);

    end
    clear volume;
    %% Crop
    cropOff = topZero+botZero;
    volume_mcorr_ChP(1:cropOff,:,:) = [];
    volume_mcorr_ChP(end-cropOff+1:end,:,:) = [];
    volume_mcorr_ChS(1:cropOff,:,:) = [];
    volume_mcorr_ChS(end-cropOff+1:end,:,:) = [];

    %% Save
    savepath = strrep(folder,'RAW DATA','Processed');
    if exist(savepath)
        savepath = savepath;
    else
        mkdir(savepath);
    end
    save(fullfile(savepath,[fn,'    _mcorr_ChP_Global']), 'volume_mcorr_ChP', '-v7.3');
    save(fullfile(savepath,[fn,'    _mcorr_ChS_Global']), 'volume_mcorr_ChS', '-v7.3');

end