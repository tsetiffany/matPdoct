%% Pigment and flow image


%% preoparing workspace and loading data
close all; clear all; clc;
% script_dir = '/ensc/IMAGEBORG/STUDENTS/DESTINY/Segmentation/Scripts/OCT_OCTA';
% addpath(script_dir);

loadloc_r = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data';
% savepath    = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data\OCTA\J001_macula_newdetection';
% savepath    = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data\OCTA\A001_9x9_WFOV';

% savepath    = 'E:\CageRPE65_ILM_LFOV';
savepath    ='E:\Cage_WT_ILM_LFOV';
% savepath    ='E:\Cage_WT_ILM_4BM';
% savepath    ='E:\JiHoon_DOPU & OCTA & avgOCT';
% savepath    = 'E:\MFOV WT NonONHRegion';
% savepath    = 'E:\Raw Data 2019\SPIE2020';
% savepath    = 'E:\CageTyr_ILM01_LFOV';

load(fullfile(loadloc_r,'cmap_dopu_r.mat')); %load the DOPU colormap
% load(fullfile(loadloc_r,'cmap_PAF_r.mat'));
load(fullfile(loadloc_r,'cmap_RPE_r.mat')); % load the RPE colourmap for en face projection of low DOPU values
load(fullfile(loadloc_r,'cmap_OCTA.mat'));
load(fullfile(loadloc_r,'cmap_OCTA_r2.mat'));


% 16_10_42-N001_3x3FOV_2   
% filename='16_53_54-CageRPE65_ILM_LFOV_0';
% filename='16_53_54-CageRPE65_ILM_LFOV_0_crop';
filename='11_14_06-Cage_WT_ILM_LFOV_crop_r';
% filename='16_57_46-CageRPE65_RPE_LFOV_2';
% filename='11_15_11-Cage_WT_RPE_LFOV';
% filename='11_15_07-Cage_WT_RPE_LFOV';
% filename='16_55_52-CageRPE65_NFL_LFOV_NEW_0';
% filename='11_14_12-Cage_WT_ILM_LFOV';

% filename='11_21_16-Cage_WT_RPE_MFOV_Upper_TiltR_adaptSmooth';
% filename='11_19_34-Cage_WT_RPE_MFOV_Upper_adaptSmooth';
% filename='11_20_21-Cage_WT_RPE_MFOV_Upper_TiltL_adaptSmooth ';

load(fullfile(savepath,[strcat(filename,'_DOPU_3_5.mat')]));
load(fullfile(savepath,[strcat(filename,'_avgOCT.mat')]));
load(fullfile(savepath,[strcat(filename,'_OCTA.mat')]));
% load(fullfile(savepath,[strcat(filename,'    _OCTA_avg.mat')]));



%% preparing OCTA & DOPU
% current process just imadjusts each B-scan, the previosu noise removal
% code was too heavy (added too much noise).  Should add histogram
% equalization across the frames as well.

% preparing OCTA
disp('Filtering OCTA...');
cplxOCTA=flipud(OCTA);
for ii=1:size(cplxOCTA,3)
    imgOCTA = mat2gray(cplxOCTA(:,:,ii));
    cplxOCTA(:,:,ii) = imadjust(imgOCTA,[0,0.2],[],2.2);
end
disp('OCTA filtered.');

% preparing DOPU
disp('Filtering DOPU...');
DOPU_filt=flipud(DOPU);

DOPU_filt(DOPU_filt>0.95) = 1; % threshold theDOPU
DOPU_test=(medfilt3(DOPU_filt, [3 5 3])); % [depth width frames]

% DOPU_test=flipud(medfilt3(DOPU));
disp('DOPU filtered.');

% size values
numPoints=size(DOPU_test,1);
numAlines=size(DOPU_test,2);
numBscans=size(DOPU_test,3);

% prepare the martices for the segmentation of the DOPU and set to zeros)
DOPU_depth_C = zeros(size(DOPU,2),size(DOPU,3));
DOPU_depth_C_bot = zeros(size(DOPU,2),size(DOPU,3));

%% Fundus image
disp('Generating fundus image...');

FUND = (imcenhance(mat2gray(squeeze(mean(avgOCT,1)))));

figure;
imagesc(FUND);colormap gray; axis equal;
title('FUND');

imwrite(FUND,fullfile(savepath,[filename,'_FUND.tif']));
% imwrite(imadjust(mat2gray(FUND_test)),fullfile(octapath,[filename,'_FUND_adapthisteq.tif']));

%% check DOPU image
botlim=5;
lowlim =10;
uplim = 20;

figure('pos',[50 50 1500 500]);
for ii=1:size(DOPU_test,3)-1 %cplxOCTA, DOPU_test
    ax1=subplot(1,2,1);
    imgA = squeeze(DOPU_test(:,:,ii));
    imagesc(imgA,[0,1]);colormap(ax1,cmap_dopu_r);colorbar;%colormap(gray, cmap_dopu);
  
    hold on;
%     plot(gcONL(:,ii),'black');
%     plot(DOPU_depth_C(:,ii)+2,'black','LineWidth', 2);
%     plot(DOPU_depth_C(:,ii)+4,'black','LineWidth', 2);
    plot(DOPU_depth_C(:,ii),'black','LineWidth', 2);
    plot(DOPU_depth_C_bot(:,ii),'black','LineWidth', 2);
    hold off;

    ax2=subplot(1,2,2);
    img2=mat2gray(flipud(squeeze(OCTA(:,:,ii))));
%     img2=(mat2gray(squeeze(cplxOCTA(:,:,ii))));
    imagesc(img2, [0, 0.1]);colormap(ax2,gray);
    title(ii);
%     
    hold on;
%     plot(gcONL(:,ii),'black');
%     plot(DOPU_depth_C(:,ii)+2,'red','LineWidth', 1);
%     plot(DOPU_depth_C(:,ii)+4,'cyan','LineWidth', 1);
    plot(DOPU_depth_C(:,ii),'red','LineWidth', 1);
    plot(DOPU_depth_C(:,ii)-botlim,'blue','LineWidth', 1);
    plot(DOPU_depth_C(:,ii)-lowlim,'yellow','LineWidth', 1);
    plot(DOPU_depth_C_bot(:,ii),'green','LineWidth', 1);
    plot(DOPU_depth_C(:,ii)-uplim,'cyan','LineWidth', 1);
    
    hold off;
    
    pause(0.0001);
%     colorbar;
end


%% Finding the cut lines (segmentation lines) of the low DOPU region
disp('Preparing method C...');

numPoints=size(DOPU_test,1);
numAlines=size(DOPU_test,2);
numBscans=size(DOPU_test,3);

startZ = 151; %depth to start searching at (usd to cut out noise)
thresh = 0.99; %threshold for finding low DOPU values (typically 0.95, lower values closer to LDR)
smoothing = 5;% if noisy, do 30, else do 5-20
botcrop = 75; % crop out any FPN near bottom

for B=1:numBscans % loop across the frames
    for A=1:numAlines-1 % loop across the A-lines
        
        % extract A-line (all depth values)
        Aline=squeeze(DOPU_test(:,A,B));
        
        % smooth out any outliers
        Aline=smoothdata(Aline,'gaussian',smoothing); 
        
        % find all indices where DOPU < thresh
        idx = find(Aline(startZ:end-botcrop)<thresh)+startZ;
        
        % extract the top and bottom values (if any are found), else set
        % the DOPU depth at the bottom of the volume
        if idx
            DOPU_depth_C(A,B)=idx(1);
            DOPU_depth_C_bot(A,B)=idx(end);
        else
            DOPU_depth_C(A,B)=numPoints;
            DOPU_depth_C_bot(A,B)=numPoints;
        end


    end
    % copy the second-last A-line just to avoid zero values
    DOPU_depth_C(end,:)=DOPU_depth_C(end-1,:);
    DOPU_depth_C_bot(end,:)=DOPU_depth_C_bot(end-1,:);
    
    % smooth out the final segmentation lines
    DOPU_depth_C(:,B)=smoothdata(DOPU_depth_C(:,B),'gaussian',50);
    DOPU_depth_C_bot(:,B)=smoothdata(DOPU_depth_C_bot(:,B),'gaussian',50);
    
    fprintf('DOPU frame: %d\n', B);
end

disp('Completed');

save(fullfile(savepath,[filename,'_DOPU_depth_C.mat']), 'DOPU_depth_C', '-v7.3');
save(fullfile(savepath,[filename,'_DOPU_depth_C_bot.mat']), 'DOPU_depth_C_bot', '-v7.3');

% displaying an elevation map

% imgDDC = imresize(notchfilter((mat2gray(DOPU_depth_C))),2);%imcenhance
imgDDC = notchfilter(imcenhance(mat2gray((DOPU_depth_C))));%imcenhance (can do log10 on the dopu depth?

% display with inverted values(smaller numbers are higher elevation originally)
figure;i1=imagesc(1-imgDDC,[0,1]);title('DOPU_depth, method C');
colormap(hot);colorbar; axis equal;

%%%%%%%%%%% check via the 'check DOPU image' section and repeat/adjust
%%%%%%%%%%% smoothing and threshold parameters until segmentation is
%%%%%%%%%%% satisfactory.


%% generate en face of low DOPU region, and also extract values like mean, thickness, an values
disp('Preparing low DOPU en face...');

% size values
numPoints=size(DOPU_test,1);
numAlines=size(DOPU_test,2);
numBscans=size(DOPU_test,3);

% parameters for the low DOPU region
lowThickness = zeros(numAlines,numBscans);
lowMean = zeros(numAlines,numBscans);
totValues=ones(size(DOPU_test));

for B=1:numBscans % loop across the frames
    for A=1:numAlines % loop across the A-lines
        
        % extract the indices (and round them)
        top = round(DOPU_depth_C(A,B));
        bot = round(DOPU_depth_C_bot(A,B));
        
        % extract the low DOPU region of inteerest
        lowROI = DOPU_test(top:bot,A,B);
        
        lowThickness(A,B) = length(lowROI); %thickness of the resion
        lowMean(A,B) = mean(lowROI); % mean value in that region
        totValues(top:bot,A,B) = lowROI; % low DOPU region values
    end
    fprintf('DOPU frame: %d\n', B);
end

lowValues = totValues(totValues~=1); % save only the low values (space independent)

% display the thickness and mean value maps
figure;
imagesc(rot90(imadjust(mat2gray(lowThickness,[0,40]),[0,1]),3),[0,1]);colormap(gray);colorbar;axis equal;
title('Thickness');
% ACTUAL max is 50, used 40 for visual clarity

figure;
imagesc(rot90(mat2gray(lowMean,[0,1]),3),[0,1]);colormap(cmap_RPE_r);colorbar;axis equal; % display wiht yellow as high dopu and blue as low
% yellow = 1-0.9, green=0.9-0.8, light blue=0.8-0.7, blue = 0.7-0.5,
% violle to magenta =0.5-0
title('DOPU value map');

% ax = gca;
% cmap_RPE_r2 = colormap(ax);
% save(fullfile(loadloc_r,'cmap_RPE_r2.mat'),'cmap_RPE_r2');

% save the parameters
save(fullfile(savepath,[filename,'_lowThickness.mat']), 'lowThickness', '-v7.3');
save(fullfile(savepath,[filename,'_lowMean.mat']), 'lowMean', '-v7.3');
save(fullfile(savepath,[filename,'_lowValues.mat']), 'lowValues', '-v7.3');


%% extract OCTA 

% determine the right segmentation area/depths using the 'cehck DOPU image'
% section (the limits are set there). set the top crop startZ to abut 40-100 or less to remove FPN

testOCTA = flipud(OCTA);
startZ=120;
% botlim = 5;
lowlim = 13;
uplim = 20;
for i = 1:size(cplxOCTA,2)
    for j = 1:size(cplxOCTA,3)
        % extract the indice (and round it)
        top = round(DOPU_depth_C(i,j));
%         top = round(gcONL(i,j));
        
%         fund_all(i,j) = mean(cplxOCTA(startZ:top-lowlim,i,j),1); 
%         fund_sup(i,j) = mean(cplxOCTA(startZ:top-uplim,i,j),1);
%         fund_deep(i,j) = mean(cplxOCTA(top-uplim:top-lowlim,i,j),1);
        
        
%         [fund_all(i,j),ia] = max(cplxOCTA(startZ:top-lowlim,i,j),[],1); 
%         [fund_sup(i,j),is] = max(cplxOCTA(startZ:top-uplim,i,j),[],1);
%         [fund_deep(i,j),id] = max(cplxOCTA(top-uplim:top-lowlim,i,j),[],1);
        
        [fund_all(i,j),ia] = max(testOCTA(startZ:top-lowlim,i,j),[],1); 
        [fund_sup(i,j),is] = max(testOCTA(startZ:top-uplim,i,j),[],1);
        [fund_deep(i,j),id] = max(testOCTA(top-uplim:top-lowlim,i,j),[],1);
%         [fund_NV(i,j),in] = max(testOCTA(top-lowlim:top-botlim,i,j),[],1);
        
        % get indices as distances from the top of the RPE
        ind_all(i,j) = top-(ia+startZ-1);
        ind_sup(i,j) = top-(is+startZ-1);
        ind_deep(i,j) = top-(id+top-uplim-1);
%         ind_NV(i,j) = top-(in+top-lowlim-1);
    end
end

figure;

OCTA_ALL = rot90(notchfilter(imcenhance(fund_all)),3);
OCTA_ALL_r = imadjust(OCTA_ALL,[0,0.5],[]);
% figure;imagesc(imadjust(OCTA_ALL));colormap(gray);axis equal;
subplot(2,2,1);imagesc(OCTA_ALL_r);colormap(gray);axis equal;
% bin = imbinarize(OCTA_ALL_r);
% figure;imagesc(bin);colormap(gray);axis equal;
title('OCTA All');

OCTA_deep = rot90(notchfilter(imcenhance(fund_deep)),3); 
OCTA_deep_r = imadjust(OCTA_deep,[0,0.5],[]);
% figure;imagesc(OCTA_deep);colormap(gray);axis equal;
subplot(2,2,2);imagesc(OCTA_deep_r);colormap(gray);axis equal;
title('OCTA deep');

OCTA_sup = rot90(notchfilter(imcenhance(fund_sup)),3); 
OCTA_sup_r = imadjust(OCTA_sup,[0,0.5],[]);
% figure;imagesc(OCTA_sup);colormap(gray);axis equal;
subplot(2,2,3);imagesc(OCTA_sup_r);colormap(gray);axis equal;
title('OCTA sup');

% combine the all (or sup) and dep for higher contrast
fund_comb = max(OCTA_ALL_r,OCTA_deep_r);
% fund_comb = OCTA_ALL_r+OCTA_deep_r;
OCTA_comb = notchfilter(imcenhance(fund_comb)); 
OCTA_comb_r = imadjust(OCTA_comb,[0,0.5],[]); 
% figure;imagesc(OCTA_comb);colormap(gray);axis equal;
subplot(2,2,4);imagesc(OCTA_comb_r);colormap(gray);axis equal;
title('OCTA comb');


% OCTA_comb_r = imadjust(OCTA_comb,[0,0.5],[],2.2); 
% % figure;imagesc(OCTA_comb);colormap(gray);axis equal;
% figure;imagesc(OCTA_comb_r);colormap(gray);axis equal;
% title('OCTA comb');

save(fullfile(savepath,[filename,'_OCTA_comb_r.mat']), 'OCTA_comb_r', '-v7.3');
imwrite(OCTA_comb_r,fullfile(savepath,[filename,'_OCTA_comb_r.tif']));

% colorMapForSaving_r = flipud(imresize(jet,[64, 3], 'nearest')); %(flip to get burgundy at top)
% figure;imagesc(ind_all,[0,50]);colormap(colorMapForSaving_r);axis equal;colorbar;
% figure;imagesc(ind_all,[0,60]);colormap(cmap_OCTA);axis equal;colorbar;
figure;imagesc(rot90(ind_all,3),[0,55]);colormap(cmap_OCTA_r2);axis equal;colorbar;


% ax = gca;
% cmap_OCTA_r2= colormap(ax);
% save(fullfile(loadloc_r,'cmap_OCTA_r2.mat'),'cmap_OCTA_r2');


save(fullfile(savepath,[filename,'_ind_all.mat']), 'ind_all', '-v7.3');


%% colour-code by depth;  ((IGNORE, EXPERIMENTING))
% grab index (or distance from the RPE) from each point where the max value
% is hit in the OCTA.

% perhaps use [M,I] = max(___), to get the value and the index I?



% combine using the lowThickness as the alpha channel;
% figure;
figure('position' , [150 150 600 500])
axBG = axes('Position',[0 0 1 1]);
blue=cat(3, zeros(size(ind_all)), zeros(size(ind_all)), ones(size(ind_all))); % overlay blue bg
b=imshow(blue); % prepare background

% red=cat(3, ones(size(ind_all)), zeros(size(ind_all)), zeros(size(ind_all))); % overlay red bg
% b=imshow(red); % prepare background

% add vessels in red
hold on;
% red=cat(3, ones(size(imgM)), zeros(size(imgM)), zeros(size(imgM))); % overlay red bg
% r=imshow(red);


axOCTA = axes('Position',[0 0 1 1]);
linkaxes([axBG,axOCTA])
% colormap(axBG , gray);

colorMapForSaving_r = (imresize(cmap_OCTA_r2,[256, 3], 'nearest'));
colormap(axOCTA , colorMapForSaving_r);

imgM = uint8(255.0 .*mat2gray(ind_all,[0,60]));
imgM_rgb = ind2rgb(imgM,colorMapForSaving_r);
m=imshow(imgM_rgb);

hold off;
mask = imresize(OCTA_comb_r,[numAlines,numBscans]); % mask intensity with OCTA
set(m, 'AlphaData', mask);

set([axBG,axOCTA],'Position',[.17 .15 .67 .815]); %from left, from botom, witdh, height
Axis3 = colorbar(axOCTA,'Position',[.1 .16 .06 .8]);
Axis3.Ticks = linspace(0, 1, 13) ; %Create 12 ticks from min to max depth
labels3 = round(linspace(0,60,13));
Axis3.TickLabels = num2cell(labels3);
Axis3.Label.String= 'Vessel distance from RPE in pixels';

%% combine en face PAF

% combine using the lowThickness as the alpha channel;
% figure;
figure('position' , [150 150 600 500])
axBG = axes('Position',[0 0 1 1]);

imshow(zeros(numAlines,numBscans)); % prepare a black background

% preparing the mean value map with colourmap of choice
colorMapForSaving = imresize(cmap_RPE_r,[256, 3], 'nearest');
imgM = rot90(uint8(255.0 .*mat2gray(lowMean,[0,1])),3);
imgM_rgb = ind2rgb(imgM,colorMapForSaving);
hold on;
axRPE=axes('Position',[0 0 1 1]);
m=imshow(imgM_rgb); % display the mean value map
% m=imagesc(mat2gray(lowMean,[0,1]),[0,1]); colormap(axRPE,cmap_RPE_r);% display the mean value map

hold off;
set(m, 'AlphaData', rot90(imadjust(mat2gray(lowThickness,[0,40]),[0,1]),3)); % mask the intensity using the thickness map
% set to 40 because the desired max of 50 is actually too dark to see
% properly.  Should realize that this is indeed 50 max instead

print(fullfile(savepath,[filename,'_PAF_DOPU.tif']),'-dtiffn') % save the tiff



% grab next axis and link together
axOCTA = axes('Position',[0 0 1 1]);
linkaxes([axBG,axRPE,axOCTA])
colormap(axRPE , cmap_RPE_r);
colormap(axBG , gray);

colorMapForSaving_r = (imresize(cmap_OCTA_r2,[256, 3], 'nearest'));
colormap(axOCTA , colorMapForSaving_r);

% add vessels from OCTA
hold on;
% red=cat(3, ones(size(imgM)), zeros(size(imgM)), zeros(size(imgM))); % overlay red bg
% r=imshow(red);
imgR = rot90(uint8(255.0 .*mat2gray(ind_all,[0,55])),3);
imgR_rgb = ind2rgb(imgR,colorMapForSaving_r);
% axOCTA=axes('Position',[0 0 1 1]);
r=imshow(imgR_rgb);
hold off;
mask = imresize(OCTA_comb_r,[numAlines,numBscans]); % mask intensity with OCTA
set(r, 'AlphaData', mask);




% show colourbars, labels, etc.
set([axBG,axRPE,axOCTA],'Position',[.17 .15 .67 .815]); %from left, from botom, witdh, height
Axis1 = colorbar(axRPE,'Position',[.1 .16 .06 .8]);
% Axis2 = colorbar(axBG,'Position',[.17 .02 .685 .0675]);
Axis2 = colorbar(axBG,'Location','southoutside','Position',[.17 .09 .67 .06]);
Axis3 = colorbar(axOCTA,'Position',[.85 .16 .06 .8]);

Axis1.Label.String= 'DOPU value';

Axis2.Ticks = linspace(0, 1, 11) ; %Create 11 ticks from min to max thickness
% labels2 = round(linspace(min(lowThickness(:)),max(lowThickness(:)),10));
labels2 = round(linspace(0,50,11));

% labels2 = round(linspace(1,54,10));


Axis2.TickLabels = num2cell(labels2);
Axis2.Label.String= 'Thickness of RPE/choroid in pixels (intensity of background)';

Axis3.Ticks = linspace(0, 1, 12) ; %Create 12 ticks from min to max depth
labels3 = round(linspace(0,55,12));
Axis3.TickLabels = num2cell(labels3);
Axis3.Label.String= 'Vessel distance from RPE in pixels';

title(filename,'Interpreter','none');

print(fullfile(savepath,[filename,'_PAF_full.tif']),'-dtiffn') % save full PAF tiff
% print(fullfile(savepath,'TEMP.tif'),'-dtiffn') % save full PAF tiff




% %re-extracting the data for MATLAB display using imagesc
% frame = getframe(gcf) ;
% save(fullfile(savepath,[filename,'_PAFframe.mat']), 'frame', '-v7.3'); % saving the frame data
% [imgD_ind, cmap]=rgb2ind(frame.cdata,256); 
% 
% % crop the image
% [row,col]=find(imgD_ind~=1);
% imgD = imgD_ind(min(row):max(row),min(col):max(col));
% 
% % % display before and after cropping
% % figure;imagesc(imgD_ind);colormap(cmap);
% % figure;imagesc(imgD);colormap(cmap);





%% histogram analysis on the en face low DOPU region

% switch based on what is desired for analysis
% strain = ', WT, volume ';
% strain = ', RPE 65';
% strain = ', Agouti';

types = {'Thickness','DOPU_values'};

for t = 1:2

%     hist_type = 'Thickness';
    % hist_type = 'DOPU_values';
    hist_type = types{t};

    if strcmp(hist_type,'Thickness')
        value = lowThickness(:);
        binlimits = [0,60];
        number = 30;
        textbox = [0.65 0.7 0.25 0.22];
    else
        value = lowValues(:);
        binlimits = [0,1];
        number = 100;
        textbox = [0.15 0.7 0.25 0.22];
    end

    % title_text = sprintf('Thickness & %s ' ,filename);
    title_text = strcat(hist_type,', volume: ',filename);


    figure;
    h_plot = histogram(value,'Normalization','probability','BinLimits', binlimits, 'NumBins',number); % set bin limits for thickness (around 80?) so the legend isn't wonky


    % Calculate the min, max, mean, median, and standard deviation
    mn=min(value);
    mx=max(value);
    me=nanmean(value);
    md=nanmedian(value);
    stdv=nanstd(value);

    % Create the labels
    minlabel=sprintf('Min - %3.2f', mn);
    maxlabel=sprintf('Max - %3.2f', mx);
    mnlabel=sprintf('Mean - %3.2f', me);
    mdlabel=sprintf('Median - %3.2f', md);
    stdlabel=sprintf('Std Deviation - %3.2f', stdv);
    % Create the textbox
    h=annotation('textbox',textbox); %distace to left, distance to top.  Adjust accordingly
    set(h,'String',{minlabel, maxlabel,mnlabel, mdlabel, stdlabel});
    title(title_text,'Interpreter', 'none');

    saveas(gcf,fullfile(savepath,[filename,'_hist_',hist_type,'.tif']));
    pause(1)
end


%% Combining histograms

% save(fullfile(savepath,'FIG4_lowValues_r.mat'), 'lowValues_r', '-v7.3');
% save(fullfile(savepath,'FIG4_lowValues_a.mat'), 'lowValues_a', '-v7.3');
% save(fullfile(savepath,'FIG4_lowThickness_r.mat'), 'lowThickness_r', '-v7.3');
% save(fullfile(savepath,'FIG4_lowThickness_a.mat'), 'lowThickness_a', '-v7.3');

% load and concatenate


hist_type = 'Thickness';
% hist_type = 'DOPU values';

if strcmp(hist_type,'Thickness')
    value_1 = lowThickness_1(:); %WT
    value_2 = lowThickness_2(:); %Agouti
    value_3 = lowThickness_3(find(lowThickness_3>1)); %albino
    binlimits = [0,50];
    number = 50;
    textbox = [0.65 0.7 0.25 0.22];
else
    value_1 = lowValues_1(:);
    value_2 = lowValues_2(:);
    value_3 = lowValues_3(:);
    binlimits = [0,1];
    number=50;
    textbox = [0.15 0.7 0.25 0.22];
end

figure;
h_plot1 = histogram(value_1,'Normalization','probability','BinLimits', binlimits, 'NumBins',number); % set bin limits for thickness (around 80?) so the legend isn't wonky
hold on
h_plot2 = histogram(value_2,'Normalization','probability','BinLimits', binlimits, 'NumBins',number);
h_plot3 = histogram(value_3,'Normalization','probability','BinLimits', binlimits, 'NumBins',number);
hold off

legend('Wild Type','Agouti','Albino','Location', 'northeast');
% legend('Wild Type','Agouti','Location', 'northeast');
% legend('Wild Type','Agouti','Albino','Location', 'northwest');
xlabel([hist_type, ' of low DOPU region']);
xlim(binlimits);
ylabel('Percent distribution');


%% Combining histograms (from MFOV nonONH region folder)

% save(fullfile(savepath,'FIG4_lowValues_r.mat'), 'lowValues_r', '-v7.3');
% save(fullfile(savepath,'FIG4_lowValues_a.mat'), 'lowValues_a', '-v7.3');
% save(fullfile(savepath,'FIG4_lowThickness_r.mat'), 'lowThickness_r', '-v7.3');
% save(fullfile(savepath,'FIG4_lowThickness_a.mat'), 'lowThickness_a', '-v7.3');

% load and concatenate
% lowThickness_r=cat(3,lowThickness_r,lowThickness);
% lowThickness_a=cat(3,lowThickness_a,lowThickness);
% lowValues_r=cat(1,lowValues_r,lowValues);
% lowValues_a=cat(1,lowValues_a,lowValues);

% hist_type = 'Thickness';
hist_type = 'DOPU values';

if strcmp(hist_type,'Thickness')
    value_rigid = lowThickness_r(:);
    value_adaptive = lowThickness_a(:);
    binlimits = [7,28];
    number = 22;
    textbox = [0.65 0.7 0.25 0.22];
else
    value_rigid = lowValues_r(:);
    value_adaptive = lowValues_a(:);
    binlimits = [0,1];
    number=22;
    textbox = [0.15 0.7 0.25 0.22];
end

figure;
h_plot1 = histogram(value_rigid,'Normalization','probability','BinLimits', binlimits, 'NumBins',number); % set bin limits for thickness (around 80?) so the legend isn't wonky
hold on
h_plot2 = histogram(value_adaptive,'Normalization','probability','BinLimits', binlimits, 'NumBins',number);
hold off

legend('Rigid','Adaptive', 'Location', 'northwest');
xlabel([hist_type, ' of low DOPU region']);
xlim(binlimits);
ylabel('Percent distribution');


%% Bar graphs, rigid vs adaptive, -3 degrees, 1 degree, 3 degrees

%%%%%%% Mean DOPU values, rigid vs adaptive, and error bounds
% degrees = categorical({'-3 degrees','1 degree','3 degrees'});
vals= [0.83 0.84; 0.82 0.83; 0.82 0.84] ;
v_errs = [0.08 0.07; 0.08 0.07; 0.08 0.07];
figure;
v = bar(vals);
xticklabels({'-3','1','3'})
hold on
% put the error bars at the center of each bar
ngroups = size(vals, 1);
nbars = size(vals, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, vals(:,i), v_errs(:,i), '.', 'Color', 'black');
end
hold off
xlabel('Directionality (degree of tilt)');
ylim([0,1.2]);
ylabel('Mean DOPU value');
legend('Rigid','Adaptive', 'Location', 'northeast');
% title('DOPU values vs directionality');
% ve = errorbar(vals,v_errs,v_errs);    
% ve.Color = [0 0 0];                            
% ve.LineStyle = 'none';  




    
%%%%%% max thickness values, rigid vs adptive
thicks = [46.00 29.00;37.00 31.00;40.00 28.00];
t_errs =[2.52 2.98;3.67 3.85;3.85 3.67];
figure;
t = bar(thicks);
xticklabels({'-3','1','3'})
hold on
% put the error bars at the center of each bar
ngroups = size(thicks, 1);
nbars = size(thicks, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, thicks(:,i), t_errs(:,i), '.', 'Color', 'black');
end
hold off
xlabel('Directionality (degree of tilt)');
ylim([0,55]);
ylabel('Max melanin region thickness');
legend('Rigid','Adaptive', 'Location', 'northeast');


%%%%%% max thickness values, rigid vs adptive (line plots)
thicks = [46.00 29.00+6.5;37.00 31.00+5.5;40.00 28.00+8];
x = [-2,0,2];
t_errs =[2.52 2.98;3.67 3.85;3.85 3.67];
figure;
tr = errorbar(x,thicks(:,1),t_errs(:,1));
% xticklabels({'-3','1','3'})
hold on
% % put the error bars at the center of each bar
% ngroups = size(thicks, 1);
% nbars = size(thicks, 2);
% % Calculating the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, thicks(:,i), t_errs(:,i), '.', 'Color', 'black');
% end
ta = errorbar(x,thicks(:,2),t_errs(:,2));
hold off
xlim([-2.5,2.5]);
xticklabels({'','-3','','','','1','','','','3',''})
xlabel('Directionality (degree of tilt)');
ylim([30,50]);
ylabel('Max melanin region thickness');
legend('Rigid','Adaptive', 'Location', 'northeast');







