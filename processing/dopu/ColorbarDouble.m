%% combine en face PAF - log scale %%
figure('position' , [150 150 600 500])
axRPE = axes;
DOPU_MAP = imagesc(20.*log10(lowMean)); colormap(cmap_RPE_Log);
imgM = uint8(255.0 .*mat2gray(lowMean,[0,1]));
set(DOPU_MAP , 'AlphaData' , imadjust(mat2gray(lowThickness)));
print(fullfile(MCHDir,[OCTAname(1 : end - 13),'_PAF_DOPU_Log.tif']),'-dtiffn') % save the tiff

axThick = axes;
linkaxes([axRPE,axThick])
colormap(axRPE , cmap_RPE_Log);
colormap(axThick , cmap_Thickness);


hold on;
red=cat(3, ones(size(imgM)), zeros(size(imgM)), zeros(size(imgM))); % overlay red bg
r=imshow(red);
hold off;
mask = imresize(OCTA_comb_r,[numAlines,numBscans]); % mask intensity with OCTA
set(r, 'AlphaData', mask);

set([axRPE,axThick],'Position',[.17 .11 .685 .815]);
Axis1 = colorbar(axRPE,'Position',[.05 .11 .0675 .815]);
Axis2 = colorbar(axThick,'Position',[.88 .11 .0675 .815]);

Axis2.Ticks = linspace(0, 1, 10) ; %Create 10 ticks from min to max thickness
labels = linspace(min(lowThickness(:)),max(lowThickness(:)),10);
Axis2.TickLabels = num2cell(labels);

print(fullfile(MCHDir,[OCTAname(1 : end - 13),'_PAF_full_Log.tif']),'-dtiffn') % save full PAF tiff

frame = getframe(gcf) ;
save(fullfile(MCHDir,[OCTAname(1 : end - 13),'_PAFframe.mat']), 'frame', '-v7.3'); % saving the frame data
[imgD_ind, cmap]=rgb2ind(frame.cdata,256); 

%%% crop the image
[row,col]=find(imgD_ind~=1);
imgD = imgD_ind(min(row):max(row),min(col):max(col));

%%% display before and after cropping
% figure;imagesc(imgD_ind);colormap(cmap);
% figure;imagesc(imgD);colormap(cmap);


%% combine en face PAF - linear scale

% combine using the lowThickness as the alpha channel;
% figure;
figure('position' , [150 150 600 500])
axBG = axes('Position',[0 0 1 1]);

imshow(zeros(numAlines,numBscans)); % prepare a black background

% preparing the mean value map with colourmap of choice
colorMapForSaving = imresize(cmap_RPE_r,[256, 3], 'nearest');
imgM = uint8(255.0 .*mat2gray(lowMean,[0,1]));
imgM_rgb = ind2rgb(imgM,colorMapForSaving);
hold on;
axRPE=axes('Position',[0 0 1 1]);
m=imshow(imgM_rgb); % display the mean value map
% m=imagesc(mat2gray(lowMean,[0,1]),[0,1]); colormap(axRPE,cmap_RPE_r);% display the mean value map

hold off;
set(m, 'AlphaData', imadjust(mat2gray(lowThickness))); % mask the intensity using the thickness map

print(fullfile(savepath,[filename,'_PAF_DOPU.tif']),'-dtiffn') % save the tiff


% grab next axis and link together
axThick = axes('Position',[0 0 1 1]);
linkaxes([axBG,axRPE,axThick])
colormap(axRPE , cmap_RPE_r);
colormap(axBG , gray);
colormap(axThick , gray);

% add vessels in red
hold on;
red=cat(3, ones(size(imgM)), zeros(size(imgM)), zeros(size(imgM))); % overlay red bg
r=imshow(red);
hold off;
mask = imresize(OCTA_comb_r,[numAlines,numBscans]); % mask intensity with OCTA
set(r, 'AlphaData', mask);


% show colourbars
set([axBG,axRPE,axThick],'Position',[.17 .11 .685 .815]);
Axis1 = colorbar(axRPE,'Position',[.05 .11 .0675 .815]);
Axis2 = colorbar(axThick,'Position',[.88 .11 .0675 .815]);

Axis2.Ticks = linspace(0, 1, 10) ; %Create 10 ticks from min to max thickness
labels = linspace(min(lowThickness(:)),max(lowThickness(:)),10);
Axis2.TickLabels = num2cell(labels);

print(fullfile(savepath,[filename,'_PAF_full.tif']),'-dtiffn') % save full PAF tiff

%re-extracting the data for MATLAB display using imagesc
frame = getframe(gcf) ;
save(fullfile(savepath,[filename,'_PAFframe.mat']), 'frame', '-v7.3'); % saving the frame data
[imgD_ind, cmap]=rgb2ind(frame.cdata,256); 

% crop the image
[row,col]=find(imgD_ind~=1);
imgD = imgD_ind(min(row):max(row),min(col):max(col));

% % display before and after cropping
% figure;imagesc(imgD_ind);colormap(cmap);
% figure;imagesc(imgD);colormap(cmap);


%%


% Original code for linking the two axis %%
ax1 = axes;
[x,y,z] = peaks;
surf(ax1,x,y,z)
view(2)
ax2 = axes;
scatter(ax2,randn(1,120),randn(1,120),50,randn(1,120),'filled')
%%Link them together
linkaxes([ax1,ax2])
%%Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'hot')
colormap(ax2,'cool')
%%Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);


%%


%% histology on the en face low DOPU region

% switch based on what is desired for analysis
title_text = sprintf('Thickness & %s ' ,upper(OCTAname(1 : end - 13)));

value = lowThickness(:);

figure;
h_plot = histogram(value,'BinLimits'  , [0 , 60] , 'Normalization','probability');

% Calculate the min, max, mean, median, and standard deviation
mn=min(value);
mx=max(value);
me=mean(value);
md=median(value);
stdv=std(value);

%%% Create legends for the data analysis
minlabel=sprintf('Min - %3.2f', mn);
maxlabel=sprintf('Max - %3.2f', mx);
mnlabel=sprintf('Mean - %3.2f', me);
mdlabel=sprintf('Median - %3.2f', md);
stdlabel=sprintf('Std Deviation - %3.2f', stdv);

% Create the textbox
h=annotation('textbox',[0.65 0.7 0.25 0.22]); %distace to left, distance to top.  Adjust accordingly
set(h,'String',{minlabel, maxlabel,mnlabel, mdlabel, stdlabel});
title(title_text , 'interpreter' , 'none' , 'FontSize' , 15);

saveas(gcf,fullfile(MCHDir,[OCTAname(1 : end - 13),'_hist','.tif']));


% end