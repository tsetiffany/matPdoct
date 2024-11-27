function corrAxialMotion(fn,volume,folder)
    %% Compute maximum cross-correlation (axial only) 
    motionA = maxxcorrAx((abs(volume).^2));
    xaxis = 1:1:size(motionA,2);

    %% Set smoothing parameter
    p = polyfit(xaxis,motionA,2);
    f = polyval(p,xaxis);

    %% Compute motion correction parameters and do motion correction       
    disp_ind = motionA - f;     
    m = size(volume,1);
    n = size(volume,3);
    topZero = max(disp_ind);
    botZero = abs(min(disp_ind));
    for k=1:n
        top = round(topZero-disp_ind(k));
        volume_mcorr(top+1:top+m,:,k) = volume(:,:,k);
    end
    clear volume;
    %% Crop
    cropOff = topZero+botZero;
    volume_mcorr(1:cropOff,:,:) = [];
    volume_mcorr(end-cropOff+1:end,:,:) = [];
    
  
    %% Save
    savepath = strrep(folder,'RAW DATA','Processed');
    if exist(savepath)
        savepath = savepath;
    else
        mkdir(savepath);
    end
    save(fullfile(savepath,[fn,'    _mcorr']), 'volume_mcorr', '-v7.3');
end