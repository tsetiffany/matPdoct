% Crop OCT volume %
clearvars

fileIdx = 1;
% Load OCT for segmentation
files_oct_A   = (dir('*-A.mat'));
files_oct_B    = (dir('*-B.mat'));
fnames_oct_A  = {files_oct_A.name}';
fnames_oct_B  = {files_oct_B.name}';
load(fnames_oct_A{fileIdx});
load(fnames_oct_B{fileIdx});
% figure(1),imagesc(imadjust(mat2gray(abs(cplxData_A(:,:,end/2))))),colormap('gray')
cplxData_A = flip(cplxData_A);
cplxData_B = flip(cplxData_B);

% crop %
depthROI        = [201, 720];
cplxData_A_crop = cplxData_A(depthROI(1):depthROI(2),:,:);
cplxData_B_crop = cplxData_B(depthROI(1):depthROI(2),:,:);

% figure(1),
% for i = 1:size(cplxData_A_crop,3)
%     imagesc(imadjust(mat2gray(abs(cplxData_A_crop(:,:,i))))),colormap('gray')
%     pause(0.1)
% end

%% save changes %%
save(fullfile(cd,fnames_oct_A{fileIdx}), 'cplxData_A_crop', '-v7.3');
save(fullfile(cd,fnames_oct_B{fileIdx}), 'cplxData_B_crop', '-v7.3');