function mat2nifti(savepath,matFile)

raw_norm = (matFile-min(matFile(:)))/(max(matFile(:)-min(matFile(:))));
raw_imadjust = (imadjustn(raw_norm));
% raw_imadjust = imadjustn(raw_norm);
raw_uint8 = uint8(255*raw_imadjust);

% save raw files% 
% disp("- saving nifiti ...")
niftiwrite(permute(raw_uint8,[2,1,3]),savepath,'Compressed',false);