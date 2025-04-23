%% Evaluating sharpness of volumes processed with adaptive kernels

close all; clear all; clc;

% loadloc = 'D:\DATA\PS-OCT';
% script_dir = 'D:\MJ\Dropbox\ProgramScripts\MatlabScripts\OCTViewer_Project';

% loadloc = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data\MS_DOPU\H001_OS_lM_PS';
% loadloc = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data\OCTA\J001_macula_newdetection';
% loadloc = 'D:\Temporary Group meeting file\A001';
% loadloc = 'E:\DNN_data\3x3mm';
% loadloc = 'E:\CageTyr_ILM01_LFOV';
% loadloc = 'E:\CageRPE65_ILM_LFOV';
% loadloc = 'E:\Cage_WT_ILM_4BM';
loadloc = 'E:\MFOV WT NonONHRegion';


script_dir = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\JonesMatrixCode_Des_Ver6_ProcForPSSAO-OCT';
addpath(script_dir);

filename = '11_21_16-Cage_WT_RPE_MFOV_Upper_TiltR_DOPU_3_5';

% non-adaptive
load(fullfile(loadloc,filename));
DOPU_filt=flipud(DOPU);
DOPU_filt(DOPU_filt>0.95) = 1; % threshold theDOPU
DOPU_test=(medfilt3(DOPU_filt, [3 5 3])); % [depth width frames]
na_DOPU = DOPU_test;

% adaptive
fn = [filename, '_adaptSmooth'];
load(fullfile(loadloc,fn));
DOPU_filt=flipud(DOPU);
DOPU_filt(DOPU_filt>0.95) = 1; % threshold theDOPU
DOPU_test=(medfilt3(DOPU_filt, [3 5 3])); % [depth width frames]
a_DOPU = DOPU_test;

%adaptive 2D (time-space)
fn = [filename, '_adaptSmooth_2D_TS'];
load(fullfile(loadloc,fn));
DOPU_filt=flipud(DOPU);
DOPU_filt(DOPU_filt>0.95) = 1; % threshold theDOPU
DOPU_test=(medfilt3(DOPU_filt, [3 5 3])); % [depth width frames]
a2_DOPU = DOPU_test;


%% loop through the frames of each and estimate the sharpness metrics
test = {na_DOPU,a_DOPU,a2_DOPU};
test_sharp = {0,0,0};
for t = 1:length(test)
    %extract the DOPU matrix
    D = test{t};
    %prepare the metrics
    BAF = [0,0,0];
    fprintf('************ DOPU test %d *************', t);
    
    for ii = 1:size(D,3)
       img =  squeeze(D(:,:,ii));
       img(isnan(img)) = 1;
       
       [B,A,F] = est_sharpness(img);
       BAF(1) = BAF(1)+B;
       BAF(2) = BAF(2)+A;
       BAF(3) = BAF(3)+F;
       fprintf('Frame processed : %d\n', ii);
    end
    % take the mean
    BAF_mean = [(BAF(1)/ii), (BAF(2)/ii), (BAF(3)/ii)];
    test_sharp{t} = BAF_mean;
end


%% display the metrics
disp(filename);
fprintf('Non-adaptive; Bisque = %f, Accutancy = %f, 20 percent high frequency = %f\n', test_sharp{1}(1),test_sharp{1}(2),test_sharp{1}(3));
fprintf('Adaptive; Bisque = %f, Accutancy = %f, 20 percent high frequency = %f\n', test_sharp{2}(1),test_sharp{2}(2),test_sharp{2}(3));
fprintf('Adaptive 2D; Bisque = %f, Accutancy = %f, 20 percent high frequency = %f\n\n', test_sharp{3}(1),test_sharp{3}(2),test_sharp{3}(3));

%% display example B-scans from the DOPUs
loadloc_r = 'C:\Users\Destiny\Documents\MATLAB\JonesMatrixCode\data';
load(fullfile(loadloc_r,'cmap_dopu_r.mat')); %load the DOPU colormap
% frame = 20;

figure;
for frame = 100
    subplot(1,3,1);
    imagesc(squeeze(na_DOPU(:,:,frame)),[0,1]);colormap(cmap_dopu_r);
    title('Non-adaptive');

    subplot(1,3,2);
    imagesc(squeeze(a_DOPU(:,:,frame)),[0,1]);colormap(cmap_dopu_r);
    title('Adaptive');

    subplot(1,3,3);
    imagesc(squeeze(a2_DOPU(:,:,frame)),[0,1]);colormap(cmap_dopu_r);
    title('Adaptive 2D');

    pause(0.05)
end