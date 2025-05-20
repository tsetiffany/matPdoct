function save_pdoct_output(avgOCT,DOPU,cmap_dopu_r,savepath,fname_save)
    disp('Saving...');
    
    patient_id = fname_save(end-5:end);
    % Frame average %
    avgNum = 5;
    for i=1:avgNum:size(avgOCT,3)
        favgOCT(:,:,ceil(i/avgNum)) = mean(avgOCT(:,:,i:i+avgNum-1),3);
        favgDOPU(:,:,ceil(i/avgNum)) = mean(DOPU(:,:,i:i+avgNum-1),3);
    end
    favgDOPU_Bscans = genDOPU_combinedBscans(favgDOPU,favgOCT,cmap_dopu_r);
    sudoEntropy = ((1-favgDOPU)); 

    %%% Save DOPU %%%
    % OCT %
    disp("- saving 3d oct ...")
    save(fullfile(savepath,[patient_id,'_raw_oct.mat']), 'avgOCT', '-v7.3');
    exportTiff(flipud(((favgOCT))), fullfile(savepath,[patient_id,'_3dv_linoct']),1)
    exportTiff(flipud(((log10(favgOCT)))), fullfile(savepath,[patient_id,'_3dv_logoct']),1)
    
    % DOPU %
    disp("- saving 3d dopu ...")
    save(fullfile(savepath,[patient_id,'_raw_dopu.mat']), 'DOPU', '-v7.3');
    exportTiff(favgDOPU_Bscans, fullfile(savepath,[patient_id,'_3dv_dopu']))
    exportTiff(flipud(sudoEntropy/0.3), fullfile(savepath,[patient_id,'_3dv_rand']))  % set color range to [0 0.3]
    
    % Save En face OCT image% 
    disp("- saving en face ...")
    f=(squeeze(mean(avgOCT,1)));
    f(f==0)=30;
    imgF=10*log10(f); 
    imgF=mat2gray(imgF);
%     imgF_test=imadjust(imgF,[0.25,0.9],[]);
%     figure(3);imshow(imgF);colormap(gray);axis equal; axis tight;
    imwrite(imgF,fullfile(savepath,[patient_id,'_enface_oct.tif']));
    close all
    
    % Save Nifti file %
    disp("- saving nifiti ...")
    mat2nifti(fullfile(savepath,[patient_id(1:5),'_oct']), flipud(favgOCT))
    mat2nifti(fullfile(savepath,[patient_id(1:5),'_dopu']), flipud(sudoEntropy))
    
%     % Save Montage_movie
%     Montage_movie_func_v2(savepath, patient_id, favgOCT,favgDOPU_Bscans,flipud(sudoEntropy/0.3));
%     saveas(gcf,fullfile(savepath,[patient_id,'_preview_dopu.tif']));
%     close all
end