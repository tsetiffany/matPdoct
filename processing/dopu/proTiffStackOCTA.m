%% tif stacking OCTA
% Requires OCTA and gcONL.  Load it before proceeding

close all; clear all; clc;
loadloc = '/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/MSM_OCTA';
script_dir = '/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/Scripts/OCT_OCTA';
addpath(script_dir);

rawpath     = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/RAW');
qcpath      = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/fundus');
mcorrpath   = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/mcorr');
octapath    = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/SV');
BVpath      = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/BVsmooth');
gcpath      = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/Segment_gc');

MSMpath      = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/MSM_OCTA');

loadloc_r = '/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/';


fn_num='16_24_25-N001_3x3FOV_3mm_MAC_G';
load(fullfile(MSMpath,[strcat(fn_num,'_OCTA.mat')]));
% load(fullfile(MSMpath,[strcat(fn_num,'_DOPU.mat')]));
% load(fullfile(MSMpath,[strcat(fn_num,'_avgOCTA.mat')]));
load(fullfile(gcpath,[fn_num,'-ONLseg.mat']));
load(fullfile(gcpath,[fn_num,'-ILMseg.mat']));
% % load(fullfile(gcpath,[fn_num,'-INLseg.mat']));

% load(fullfile(loadloc_r,'cmap_dopu_r.mat')); %load the DOPU colormap
% load(fullfile(loadloc_r,'cmap_PAF_r.mat'));

%% Tif-stacking OCTA
% prepare the OCTA
cplxOCTA=flipud(OCTA);
for ii=1:size(cplxOCTA,3)
    imgOCTA = mat2gray(cplxOCTA(:,:,ii));
    imgOCTA = imadjust(imgOCTA,[0, 0.1],[]);
    for i = 1:size(imgOCTA,2)
        nOCT(i) = max(imgOCTA(1:100,i)); % orig 1:20
        nvOCT(i) = var(imgOCTA(1:100,i)); % orig 1:20
    end
    OCT_N = median(nOCT);
    OCT_Nvar=median(nvOCT);
%     cplxOCTA(:,:,ii) = imadjust(imgOCTA-OCT_N, [0, 0.1],[],2.9);
%     cplxOCTA(:,:,ii) = wiener2(cplxOCTA(:,:,ii), [3, 3],OCT_Nvar);
    cplxOCTA(:,:,ii) = wiener2(imgOCTA, [5, 5],OCT_Nvar);
    cplxOCTA(:,:,ii) = imadjust(cplxOCTA(:,:,ii)-OCT_N, [0, 0.3],[],2.0);

    %     imgOCTA=medfilt2(imgOCTA);
    fprintf('Frame: %d\n', ii);
end

%% Check the ONL
figure;
for ii = 1:500 %iterate over # of frames
%     imagesc(imadjust(mat2gray(flipud(((avgOCT(:,:,ii)))),[0,0.5])));
%     imagesc(flipud(mat2gray(avgOCT(:,:,ii))),[0 0.2]);
    imagesc((mat2gray(cplxOCTA(:,:,ii))),[0 0.5]);
    colormap(gray);

    hold on %plot the lines over the desired image
    plot(gcILM(:,ii)+5, 'green');
%     plot(gcINL(:, ii), 'red');
    plot(gcONL(:, ii)-15, 'cyan');
%     plot(gcONL(:, ii)-75, 'cyan');
    
    hold off
    title(fn_num);
    pause(0.01)
end

%% Cut out the RNFL and below the ONL
% Zero everything above the RNFL and below the ONL

for i=1:size(cplxOCTA,3) %loop across the frames
    for j=1:size(cplxOCTA,2) % loop across the A-lines
       cplxOCTA(1:gcILM(i,j)+5,i,j)=0;
       cplxOCTA(gcONL(i,j)-15:end,i,j)=0;
    end
end

%% Write the tif to a tif file
% savepath=octapath;
savepath='C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data\OCTA\H001_PAF_test';
outputFileName = fullfile(savepath,[fn_num,'_OCTA.tif']);


for K=80:160
%    imgOCTA = notchfilter(imcenhance(squeeze(cplxOCTA(K, :, :))));
   imgOCTA = notchfilter(squeeze(cplxOCTA(K, :, :)));
   imwrite(imgOCTA, outputFileName, 'WriteMode', 'append','Compression','none');
   fprintf('Frame: %d\n', K);
end


