clearvars
addpath(genpath('/Users/YM/Documents/GitHub/COIL/Project-PD_OCT/Melanoma/code'));
cd('C:\Users\coil_\Desktop\Github\Project-PD_OCT\Melanoma\code\JonesMatrixCode\JonesMatrixCode_Des_Ver6_ProcForPSSAO-OCT')
%%
fileloc = 'C:\Users\coil_\OneDrive - UBC';
segloc = 'C:\Users\coil_\OneDrive - UBC\Share\PDOCT-Segmentation';
filename = 'PS023_OS_OCT';

vol = load(fullfile(fileloc,'11_17_29-_DOPU_Bscans_3_5_dopu.mat')).DOPU_test;
seg = double(permute(niftiread(fullfile(segloc,'(post interpol) PS023_OS_seg.nii.gz')),[2,1,3]));
entropy = 1 - flipud(vol);

frame = 200;
figure,imshow(entropy(:,:,frame),[])
figure,imshow(entropy(:,:,frame).*seg(31:530,:,frame),[])

rDOPU = entropy.*seg(31:530,:,:);
cDOPU = entropy.* (1-seg(31:530,:,:));

rDOPU_sm = imgaussfilt(rDOPU,2);
cDOPU_sm = imgaussfilt(cDOPU,2);

figure(1),imshow((squeeze(mean(rDOPU,1))),[]),colormap('hot')
figure(2),imshow(squeeze(mean(cDOPU,1)),[]),colormap('hot')

rDOPU = zeros(size(DOPU_depth_C));
rDOPU_vol = entropy;
for i = 1:1000
    for j = 1:1000
        upperlim = round(DOPU_depth_C_sm(i,j));
        lowerlim = upperlim + rHite;
        lowerlim(lowerlim>size(entropy,1)) = 1;
        upperlim(lowerlim>size(entropy,1)) = 1;
        rDOPU(i,j) = 1-mean(entropy,1);
%         rDOPU(i,j) = 1-mean(avgOCT(upperlim:lowerlim,i,j),1);  
        rDOPU_vol(1:upperlim-1,i,j) = 0;
        rDOPU_vol(lowerlim+1:end,i,j) = 0;
    end
    if mod(i,100) == 0
        i
    end
end
figure,imshow(1-rDOPU,[]),colormap('hot')

savepath = '/Users/YM/Desktop/';
exportTiff(rDOPU_vol, fullfile(savepath,['rDOPU']))  