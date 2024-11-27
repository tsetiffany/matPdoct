%% Pigment and flow image
% Requires DOPU, OCTA and gcILM.  Load them before proceeding

close all; clear all; clc;
% loadloc = '/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/MSM_OCTA';
% script_dir = '/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/Scripts/OCT_OCTA';
% addpath(script_dir);

% rawpath     = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/RAW');
% qcpath      = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/fundus');
% mcorrpath   = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/mcorr');
% octapath    = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/SV');
% BVpath      = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/BVsmooth');
% gcpath      = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/Segment_gc');
% 
% MSMpath      = fullfile('/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/MSM_OCTA');

% loadloc_r = '/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/DATA/';
loadloc_r = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data';
% savepath    = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data\OCTA\J001_macula_newdetection';
% savepath    = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data\OCTA\A001_9x9_WFOV';
savepath    = 'E:\CageRPE65_ILM_LFOV';
% savepath    ='E:\Cage_WT_ILM_LFOV';



load(fullfile(loadloc_r,'cmap_dopu_r.mat')); %load the DOPU colormap
load(fullfile(loadloc_r,'cmap_PAF_r.mat'));
% 15_05_48-N001_3x3FOV_3mm_MAC
% 15_01_43-H001_3x3FOV_3mm_MAC_r
% 15_04_25-H001_6x6FOV_2mm_MAC
% 15_05_48-N001_3x3FOV_3mm_MAC_r (DOPU kernel 1,5)
% 16_00_48-N001_6x6FOV_2mm_MAC

% 16_10_42-N001_3x3FOV_2   
filename='16_53_54-CageRPE65_ILM_LFOV_0';
load(fullfile(savepath,[strcat(filename,'_OCTA.mat')]));
load(fullfile(savepath,[strcat(filename,'_DOPU_3_5.mat')]));
load(fullfile(savepath,[strcat(filename,'_avgOCT.mat')]));
% load(fullfile(savepath,[strcat(filename,'_OCTA_avg.mat')]));
% load(fullfile(gcpath,[filename,'-ONLseg.mat']));
% load(fullfile(gcpath,[filename,'-ILMseg.mat']));
% load(fullfile(gcpath,[filename,'-INLseg.mat']));


%% preparing OCTA
disp('Filtering OCTA...');
cplxOCTA=flipud(OCTA);
for ii=1:size(cplxOCTA,3)
% %     cplxOCTA(:,:,ii)=imadjust(mat2gray(squeeze(cplxOCTA(:,:,ii))),[0, 0.05], []);
% %     cplxOCTA(:,:,ii)= adapthisteq(mat2gray(cplxOCTA(:,:,ii)),'NBins',200,'ClipLimit',0.009);
% %     cplxOCTA(:,:,ii)= adapthisteq(mat2gray(cplxOCTA(:,:,ii)));
%     %%%%%%%%%%%%% testing out noise cancellation
%     imgOCTA = mat2gray(cplxOCTA(:,:,ii));
%     imgOCTA = imadjust(imgOCTA,[0, 0.1],[]);
%     for i = 1:size(imgOCTA,2)
%         nOCT(i) = max(imgOCTA(1:20,i)); % orig 1:20
%         nvOCT(i) = var(imgOCTA(1:20,i)); % orig 1:20
%     end
%     OCT_N = median(nOCT);
%     OCT_Nvar=median(nvOCT);
% %     cplxOCTA(:,:,ii) = imadjust(imgOCTA-OCT_N, [0, 0.1],[],2.9);
% %     cplxOCTA(:,:,ii) = wiener2(cplxOCTA(:,:,ii), [3, 3],OCT_Nvar);
%     cplxOCTA(:,:,ii) = wiener2(imgOCTA, [3, 3],OCT_Nvar);
%     cplxOCTA(:,:,ii) = imadjust(cplxOCTA(:,:,ii)-OCT_N, [0, 0.1],[],2.9);
% %     imgOCTA=medfilt2(imgOCTA);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    imgOCTA = mat2gray(cplxOCTA(:,:,ii));
    cplxOCTA(:,:,ii) = imadjust(imgOCTA,[0,0.2],[],2.2);
end
disp('OCTA filtered.');


%% preparing DOPU
disp('Filtering DOPU...');
DOPU_filt=flipud(DOPU);
% threshold theDOPU
DOPU_filt(DOPU_filt>0.95) = 1;
DOPU_test=(medfilt3(DOPU_filt, [3 5 3]));

% DOPU_test=flipud(medfilt3(DOPU));
disp('DOPU filtered.');
DOPU_depth_C = zeros(size(DOPU,2),size(DOPU,3));
DOPU_depth_C_bot = zeros(size(DOPU,2),size(DOPU,3));
%%
loadloc_r = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data';
load(fullfile(loadloc_r,'cmap_dopu_r.mat')); %load the DOPU colormap

%% check DOPU image
figure('pos',[50 50 1000 500]);
for ii=1:size(DOPU_test,3)-1 %cplxOCTA, DOPU_test
    ax1=subplot(1,2,1);
    imgA = squeeze(DOPU_test(:,:,ii));
    imagesc(imgA);colormap(ax1,cmap_dopu_r);%colormap(gray, cmap_dopu);
  
    hold on;
%     plot(gcONL(:,ii),'black');
%     plot(DOPU_depth_C(:,ii)+2,'black','LineWidth', 2);
%     plot(DOPU_depth_C(:,ii)+4,'black','LineWidth', 2);
    plot(DOPU_depth_C(:,ii),'black','LineWidth', 2);
    plot(DOPU_depth_C_bot(:,ii),'black','LineWidth', 2);
    hold off;

    ax2=subplot(1,2,2);
%     img2=flipud(squeeze(avgOCT(:,:,ii)));
    img2=mat2gray(squeeze(cplxOCTA(:,:,ii)));
    imagesc(img2, [0, 0.5]);colormap(ax2,gray);
    title(ii);
    
    hold on;
%     plot(gcONL(:,ii),'black');
%     plot(DOPU_depth_C(:,ii)+2,'red','LineWidth', 1);
%     plot(DOPU_depth_C(:,ii)+4,'cyan','LineWidth', 1);
    plot(DOPU_depth_C(:,ii)-10,'red','LineWidth', 1);
    plot(DOPU_depth_C_bot(:,ii),'cyan','LineWidth', 1);
    plot(DOPU_depth_C(:,ii)-30,'green','LineWidth', 1);
    hold off;
    
    pause(0.01);
%     colorbar;
end

%% Fundus image
disp('Generating fundus image...');

FUND = imcenhance(mat2gray(squeeze(mean(avgOCT,1))));

figure;
imagesc(FUND);colormap gray; axis equal;
title('FUND');

imwrite(FUND,fullfile(savepath,[filename,'_FUND.tif']));
% imwrite(imadjust(mat2gray(FUND_test)),fullfile(octapath,[filename,'_FUND_adapthisteq.tif']));

%% Check OCTA
figure;
for ii = 1:size(cplxOCTA,3) %iterate over # of frames
    imagesc(imadjust((mat2gray(filtavgOCT(:,:,ii)))));
    colormap(gray);
    
    hold on %plot the lines over the desired image
    plot(gcILM(:,ii)+5, 'green');
    plot(gcINL(:, ii), 'red');
    plot(gcONL(:, ii)+18, 'cyan');
    plot(gcONL(:, ii)+48, 'cyan');
    
    hold off
    title(ii);
    pause(0.01)
end

%% OCTA en-face
for i = 1:size(cplxOCTA,2)
    for j = 1:size(cplxOCTA,3)
        fund_all(i,j) = mean(cplxOCTA(gcILM(i,j)+10:gcONL(i,j)-15,i,j),1);
        fund_deep(i,j) = mean(cplxOCTA(gcINL(i,j):gcONL(i,j)-15,i,j),1); %ONL was  - 20
        fund_sup(i,j) = mean(cplxOCTA(gcILM(i,j)+5:gcINL(i,j),i,j),1); %INL was +10
        
        fund_chor(i,j) = mean(filtavgOCT(gcONL(i,j)+18:gcONL(i,j)+48,i,j),1);

    end
end
% OCTA_ALL = notchfilter(imcenhance(fund_all));
% figure;imagesc(OCTA_ALL);colormap(gray);axis equal;
% OCTA_deep = notchfilter(imcenhance(fund_deep)); 
% figure;imagesc(OCTA_deep);colormap(gray);axis equal;
% OCTA_sup = notchfilter(imcenhance(fund_sup)); 
% figure;imagesc(OCTA_sup);colormap(gray);axis equal;
OCTA_chor = notchfilter(imcenhance(fund_chor)); 
figure;imagesc(OCTA_chor);colormap(gray);axis equal;




%% Method A) alpha filtering of OCTA on top of the depth map
% % May need to adjust the colormap or try and plot with two colormaps
% 
% imgDD = notchfilter(imcenhance(mat2gray(DOPU_depth)));
% figure;i1=imagesc(imadjust(imgDD));title('DOPU_depth');
% colormap(colorcube);colorbar;
% hold on;
% i2=imagesc(imadjust((OCTA_ALL)));%colormap(i2,gray);colorbar;
% alpha(i2,'color');
% hold off
% 
%% Method B) surface plot of depth map over a 2D plot of the OCTA
% % need to turn off edge color and maybe smooth out the values of the depth
% % maybe try median filtering the depth values, or smooth it B-scan by
% % B-scan
% figure;
% surf(imgDD,'EdgeColor','none');
% hold on;
% i2=imagesc(imadjust(OCTA_ALL));%colormap(i2,gray);
% hold off;




%% Method C), don't segment from top of retina, do from top of volume.  Take
% median value where DOPU <1 as the depth.
disp('Preparing method C...');

numPoints=size(DOPU_test,1);
numAlines=size(DOPU_test,2);
numBscans=size(DOPU_test,3);

startZ = 100; %depth to start searching at
thresh = 0.95; %threshold for finding low DOPU values

for B=1:numBscans % loop across the frames
    for A=1:numAlines % loop across the A-lines
        % extract A-line (all depth values)
        Aline=squeeze(DOPU_test(:,A,B));
        
        % remove outliers along the A-line by a moving mean window
%         Aline(isoutlier(Aline,'movmean',50))=1;
%         Aline=smooth(1:length(Aline),Aline,0.1,'rloess'); %smooth, robust lossless, 10% span
        Aline=smoothdata(Aline,'gaussian',10); % if noisy, do 30, else do 10
        
        % find all indices where DOPU < thresh
        idx = find(Aline(startZ:end-0)<thresh)+startZ;
%         % get the median value
%         DOPU_depth_C(A,B)=median(idx(1:end));

        if idx
            DOPU_depth_C(A,B)=idx(1);
            DOPU_depth_C_bot(A,B)=idx(end);
        else
            DOPU_depth_C(A,B)=numPoints;
            DOPU_depth_C_bot(A,B)=numPoints;
        end
%         DOPU_depth_C(A,B)=median(find(squeeze(DOPU_test(120:end,A,B))<1));

    end
    % remove outliers along the A-lines by a moving mean window
%     DOPU_depth_C(isoutlier(DOPU_depth_C(:,B),'movmean',100),B)=NaN;
    DOPU_depth_C(:,B)=smoothdata(DOPU_depth_C(:,B),'gaussian',50);
    DOPU_depth_C_bot(:,B)=smoothdata(DOPU_depth_C_bot(:,B),'gaussian',50);
    
    fprintf('DOPU frame: %d\n', B);
end
disp('Completed');

save(fullfile(savepath,[filename,'_DOPU_depth_C.mat']), 'DOPU_depth_C', '-v7.3');
save(fullfile(savepath,[filename,'_DOPU_depth_C_bot.mat']), 'DOPU_depth_C_bot', '-v7.3');

% % Get rid of NANs
% DOPU_depth_C(isnan(DOPU_depth_C))=0; % setting to zero
% DOPU_depth_Cr=DOPU_depth_C;

% % Interpolating the NANs
% nanx = isnan(DOPU_depth_C);
% t    = 1:numel(DOPU_depth_C);
% DOPU_depth_C(nanx) = interp1(t(~nanx), DOPU_depth_C(~nanx), t(nanx));
% % Catch NaNs at endpoints
% nane = find(isnan(DOPU_depth_C)==1);
% if length(nane)
%     for jj = 1:length(nane)
%         if nane(jj)==1
%             DOPU_depth_C(nane(jj)) = DOPU_depth_C(nane(jj)+1);
%         else
%             DOPU_depth_C(nane(jj)) = DOPU_depth_C(nane(jj)-1);
%         end
%     end
% end

% imgDDC = imresize(notchfilter((mat2gray(DOPU_depth_C))),2);%imcenhance
imgDDC = notchfilter(imcenhance(mat2gray((DOPU_depth_C))));%imcenhance (can do log10 on the dopu depth?

% display with inverted (smaller numbers are higher elevation)
figure;i1=imagesc(1-imgDDC,[0,1]);title('DOPU_depth, method C');
colormap(hot);colorbar; axis equal;
% hold on;
% i2=imagesc(imadjust((OCTA_ALL)));%colormap(i2,gray);colorbar;
% alpha(i2,'color');axis equal;
% hold off

% OCTA_ALL_bin=imbinarize(OCTA_ALL,'adaptive');
% figure;imshow(OCTA_ALL_bin);

%% generate en face of low DOPU region, and also extract values like mean, thickness, an values
disp('Preparing low DOPU en face...');

numPoints=size(DOPU_test,1);
numAlines=size(DOPU_test,2);
numBscans=size(DOPU_test,3);

lowThickness = zeros(numAlines,numBscans);
lowMean = zeros(numAlines,numBscans);
totValues=ones(size(DOPU_test));

for B=1:numBscans % loop across the frames
    for A=1:numAlines % loop across the A-lines
        
        % extract the indices (and round them)
        top = round(DOPU_depth_C(A,B));
        bot = round(DOPU_depth_C_bot(A,B));
        
        % extract the low DOPU regoin of inteerest
        lowROI = DOPU_test(top:bot,A,B);
        
        lowThickness(A,B) = length(lowROI); %thickness of the resion
        lowMean(A,B) = mean(lowROI); % mean value in that region
%         totValues = [totValues, lowROI']; % save these values (space independent)
        totValues(top:bot,A,B) = lowROI;
    end

    fprintf('DOPU frame: %d\n', B);
end

lowValues = totValues(totValues~=1);

figure;
imagesc(imadjust(mat2gray(lowThickness)));colormap gray;colorbar;

figure;
imagesc(mat2gray(1-lowMean));colormap(parula);colorbar; % display wiht yellow as high dopu and blue as low


save(fullfile(savepath,[filename,'_lowThickness.mat']), 'lowThickness', '-v7.3');
save(fullfile(savepath,[filename,'_lowMean.mat']), 'lowMean', '-v7.3');
save(fullfile(savepath,[filename,'_lowValues.mat']), 'lowValues', '-v7.3');


%% extract OCTA 
for i = 1:size(cplxOCTA,2)
    for j = 1:size(cplxOCTA,3)
        % extract the indice (and round them)
        top = round(DOPU_depth_C(A,B));
        
        fund_all(i,j) = mean(cplxOCTA(1:top-10,i,j),1);
        
        fund_sup(i,j) = mean(cplxOCTA(1:top-30,i,j),1);
        fund_deep(i,j) = mean(cplxOCTA(top-30:top-10,i,j),1);
        
        
    end
end
OCTA_ALL = notchfilter(imcenhance(fund_all));
OCTA_ALL_r = imadjust(OCTA_ALL,[0,0.7],[]);
% figure;imagesc(imadjust(OCTA_ALL));colormap(gray);axis equal;
figure;imagesc(OCTA_ALL_r);colormap(gray);axis equal;
% bin = imbinarize(OCTA_ALL_r);
% figure;imagesc(bin);colormap(gray);axis equal;

OCTA_deep = notchfilter(imcenhance(fund_deep)); 
OCTA_deep_r = imadjust(OCTA_deep,[0,0.7],[]);
% figure;imagesc(OCTA_deep);colormap(gray);axis equal;
figure;imagesc(OCTA_deep_r);colormap(gray);axis equal;

% OCTA_sup = notchfilter(imcenhance(fund_sup)); 
% figure;imagesc(OCTA_sup);colormap(gray);axis equal;

fund_comb = max(OCTA_ALL_r,OCTA_deep_r);
% fund_comb = fund_sup+fund_deep;
OCTA_comb = notchfilter(imcenhance(fund_comb)); 
OCTA_comb_r = imadjust(OCTA_comb,[0,0.7]); 
% figure;imagesc(OCTA_comb);colormap(gray);axis equal;
figure;imagesc(OCTA_comb_r);colormap(gray);axis equal;

%% combines en face PAF

% combine using the lowThickness as the alpha channel;
figure;
imshow(zeros(numAlines,numBscans));

% colorMapForSaving = imresize(cmap_dopu_r,[256, 3], 'nearest');
% imgM = uint8(255.0 *mat2gray(lowMean)); %without mat2gray the images lose resolution and look like gray2ind, even though mat2gray on ingD doesn't seem to change teh values
% imgM_rgb=ind2rgb(floor(imgM),colorMapForSaving);

colorMapForSaving = imresize(parula,[256, 3], 'nearest');
imgM = uint8(255.0 *mat2gray(1-lowMean));
imgM_rgb = ind2rgb(imgM,colorMapForSaving);
hold on;
m=imshow(imgM_rgb);
hold off;
set(m, 'AlphaData', imadjust(mat2gray(lowThickness)));

print(fullfile(savepath,[filename,'_PAF_DOPU.tif']),'-dtiffn')

% add vessels in red
hold on;
red=cat(3, ones(size(imgM)), zeros(size(imgM)), zeros(size(imgM)));
r=imshow(red);
hold off;
mask = imresize(OCTA_comb_r,[numAlines,numBscans]);
set(r, 'AlphaData', mask);

print(fullfile(savepath,[filename,'_PAF_full.tif']),'-dtiffn')

frame = getframe(gcf) ;
save(fullfile(savepath,[filename,'_PAFframe.mat']), 'frame', '-v7.3');

[imgD_ind, cmap]=rgb2ind(frame.cdata,256); % 64 is the highest number allowable if the bg is black

% crop the image
[row,col]=find(imgD_ind~=1);


imgD = imgD_ind(min(row):max(row),min(col):max(col));

figure;imagesc(imgD_ind);colormap(cmap);

figure;imagesc(imgD);colormap(cmap);

% imshow(frame.cdata);






%% histology on the en face low DOPU region
value = lowThickness(:);
% value = lowValues(:);
title_text = 'Thickness, WT';
% title_text = 'DOPU values, WT';

h_plot = histogram(value,'Normalization','probability');

% Calculate the min, max, mean, median, and standard deviation
mn=min(value);
mx=max(value);
me=mean(value);
md=median(value);
stdv=std(value);
% Create the labels
minlabel=sprintf('Min - %3.2f', mn);
maxlabel=sprintf('Max - %3.2f', mx);
mnlabel=sprintf('Mean - %3.2f', me);
mdlabel=sprintf('Median - %3.2f', md);
stdlabel=sprintf('Std Deviation - %3.2f', stdv);
% Create the textbox
h=annotation('textbox',[0.7 0.8 0.1 0.1]); %distace to left, distance to top, 
set(h,'String',{minlabel, maxlabel,mnlabel, mdlabel, stdlabel});
title(title_text);

%% Method D: weighted linearly along depth;
disp('Preparing method D...');

dd=1:size(DOPU,1); %linear depth weights
% dd=size(DOPU,1):-1:1; %linear depth weights in reverse order
dd_sum=sum(dd);% default total sum value

for B=1:size(DOPU,3) % loop across the frames
    for A=1:size(DOPU,2) % loop across the A-lines
        
        % multiply each DOPU value against equivalent depth
        DOPU_depth_D(A,B)=sum(squeeze(DOPU_test(:,A,B)).*dd')./dd_sum;

    end    
    fprintf('DOPU frame: %d\n', B);
end
disp('Completed');
% imgDDD = imresize(notchfilter((mat2gray(DOPU_depth_D))),2);%imcenhance
imgDDD = (imcenhance(mat2gray(DOPU_depth_D)));%imcenhance
% imgDDD

figure;i1=imagesc(1-imgDDD);title('DOPU_depth, method D');
colormap(cmap_dopu_r);colorbar; axis equal;
% hold on;
% i2=imagesc(imadjust((OCTA_ALL)));%colormap(i2,gray);colorbar;
% alpha(i2,'color');axis equal;
% hold off



%% Log scale adjustment

% filter the avgOCT
disp('Filtering avgOCT...');
filtavgOCT=mat2gray(flipud(avgOCT));
for ii=1:size(filtavgOCT,3)
    %%%%%%%%%%%%% testing out noise cancellation
    imgavgOCT = (filtavgOCT(:,:,ii));
%     imgavgOCT = imadjust(imgavgOCT,[0, 0.1],[]);
    for i = 1:size(imgavgOCT,2)
        nOCT(i) = max(imgavgOCT(1:20,i)); % orig 1:20
        nvOCT(i) = var(imgavgOCT(1:20,i)); % orig 1:20
    end
    OCT_N = median(nOCT);
    OCT_Nvar=median(nvOCT);
    filtavgOCT(:,:,ii) = wiener2(imgavgOCT, [3, 3],OCT_Nvar);
%     filtavgOCT(:,:,ii) = imadjust(cplxavgOCT(:,:,ii)-OCT_N, [0, 0.1],[],2.9);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% B-scan frame
ii=220;

% take the logscale of the avgOCT
imgL=20.*log10(filtavgOCT(:,:,ii));
imgL=mat2gray(imgL);

imgL_test=imadjust(imgL,[0.2,0.7],[],1.1);
% imgL_test=imadjust(imgL,[0.3,0.8],[],1.9);

% figure;subplot(1,2,1);imagesc(imgL);colormap(gray);axis equal;
% subplot(1,2,2);
% imgL_test=imadjust(mat2gray(flipud(OCTA(:,:,250))));
figure;imagesc(imgL_test);colormap(gray);axis equal;

% convert the DOPU image to RGB then HSV
imgD=(DOPU_test(:,:,ii));
figure;imagesc(imgD);colormap(cmap_dopu_r);axis equal;

% colorMapForSaving = imresize(cmap_dopu_r,[65536, 3], 'nearest');
% imgD_rgb=ind2rgb(ceil(imgD.*65535),colorMapForSaving);
% imgD_hsv=rgb2hsv(imgD_rgb); %h=1, s=2, v=3
% 
% % replace the V channel with the log image
% imgD_hsv(:,:,3)=imgL_test;
% imgD_rgb_r=hsv2rgb(imgD_hsv);
% 
% % display the images
% figure;subplot(2,1,1);imshow(imgD_rgb);
% subplot(2,1,2);imshow(imgD_rgb_r);



% %% preparing PAF
% disp('Preparing PAF...');
% for B=1:size(DOPU,3) % loop across the frames
%     for A=1:size(DOPU,2) % loop across the A-lines
%         surface=gcILM(A,B)+7; % surface index
%         idx=1;
%         D_sum=ones(size(DOPU,1)-surface+1,1);
%         for z=surface:size(DOPU,1) % loop down the depth starting from top of retina
%             D_sum(idx)=exp(-(z-surface).*sum(1-DOPU_test(surface:z,A,B)));
%             if D_sum(idx)<0.001
%                 break;
%             end
%             idx=idx+1;
%         end
%         paf(A,B)=sum(squeeze(cplxOCTA(surface:size(DOPU,1),A,B)).*(D_sum(:).^0.01));
%         
% %         DOPU_depth(A,B)=sum(D_sum(:).^0.01);
%         DOPU_depth(A,B)=idx;
%     end    
%     
%     fprintf('PAF frame: %d\n', B);
% end
% disp('PAF completed, generating en face');
% % save(fullfile(loadloc,[fn_num,'_paf.mat']), 'paf', '-v7.3');
% %% en-face PAF
% imgDD = notchfilter(imcenhance(1-mat2gray(DOPU_depth)));
% % imgPAF=imadjust(imgPAF,[],[],2.2);
% figure;imagesc(imadjust(imgDD));title('DOPU_depth');colormap(jet);colorbar;
% 
% % imwrite(imadjust(mat2gray(imgPAF)),fullfile(loadloc,['\',fn_num,'_PAF.tif']));

%% Choroidal flow
% performed only on wide field 9x9 images.  Requires: RPE segmented (eg.
% from the elevation map for PAF) and avgOCT data.

% take the en face of volume at each depth from RPE down to choroid bottom

% do adaptive thresholding to segment out the vessels

% combine the layers together to form a depth-resolved vessel map

% how to overlay on RPE elevation???  Some sort of blending?  
choroid = zeros(60,size(avgOCT,2),size(avgOCT,3));

for ii=1:size(avgOCT,3) %loop throguh frames
    for jj=1:size(avgOCT,2) %loop across B-scans
        for kk=1:60 %loop 60 pixels down
            depth = round(DOPU_depth_C(jj,ii))+(kk-1);
            if depth > size(avgOCT,1)
                depth = size(avgOCT,1);
            end
            choroid(kk,jj,ii) = avgOCT(depth,jj,ii);
        end
    end
end

ef_chr = mat2gray(squeeze(mean(choroid,1)));
figure;imagesc(imadjust(ef_chr));axis equal; title('En-face choroid');colormap(gray);colorbar;
