%% OCT/OCTA Processing Pipeline %%
% updated : 2023.04.08
% tic
% clearvars, clc, close all
%function OCT_Process_VCSEL_Sockeye(folder,savepath,fileIdx,dispROI)
%%
%%-----------------------------------------------------------------------------------------------------%%
% local computer path % %
tic
datapath = 'I:\15May25_MJ_OD';
code_loc = 'C:\Users\tiffa\Documents\1. Projects\PDOCT Processing\main_vctrl\matPdoct_v1.0';
addpath(genpath(code_loc));
fileIdx = 1;
savepath = datapath;
dispROI = [30, 500];

% cluster path %
% SCRIPTLOC = '/arc/home/gracesoo/Code/matPdoct';
% addpath(genpath(SCRIPTLOC));
%%-----------------------------------------------------------------------------------------------------%%

dispMaxOrder    = 5;
% fileIdx = 1;

%%%%% Find filenames of RAW files to be processed %%%%%
cd(datapath);
files   = (dir('*.unp'));
fnames  = {files.name}';
if fileIdx > length(fnames)
    error('Job terminated.\n Index %d exeeds maximum number of files.', fileIdx)
end
fn = fnames{fileIdx};

%%% Load Acquisition Parameters %%%
parameters    = getParameters(fn);
numPoints     = parameters(1)*2;
numAscans     = parameters(2)/2;
numBscans     = parameters(3);
numCscans     = parameters(4);
numMscans     = parameters(5);

%%% LUT %%%
fn_ResParam = 'LUT_5050.bin';
fid_ResParam = fopen(fn_ResParam);
rescaleParam = fread(fid_ResParam, 'double');
LUT =  rescaleParam;
fclose(fid_ResParam);

% Load colormap %
cmap_dopu = load('cmap_dopu_r.mat').cmap_dopu_r;    %load the DOPU colormap

% Create ProcdData\PS\OD or OS folders %
cd(savepath)
[~,fname_save,~] = fileparts(fn);
pathSplit = regexp(datapath,'\','split');
if length(pathSplit)==1
    pathSplit = regexp(datapath,'/','split');
end
% patientId = pathSplit{end-1};
% side = pathSplit{end-1};
if ~isfolder('ProcdData')
    mkdir ProcdData
end
cd('ProcdData') %PS...
if ~isfolder(pathSplit{end-1})
    mkdir(pathSplit{end-1})
end
cd(pathSplit{end-1}) %OD OS
if ~isfolder(pathSplit{end})
    mkdir(pathSplit{end})
end
cd(pathSplit{end}) 
if ~isfolder(fname_save)
    mkdir(fname_save)
end
cd(fname_save)
if ~isfolder('log')
    mkdir('log')
end
cd(datapath)
process_path = fullfile(savepath,'ProcdData',pathSplit{end-1},pathSplit{end},fname_save);
log_path = fullfile(savepath,'ProcdData',pathSplit{end-1},pathSplit{end},fname_save,'log');
file_id      = char(fname_save(end-5:end));

ref_frame_disp = round(numBscans/2);
ref_frame_noise = 2;%numBscans - 1;
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPERSION ESTIMATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Dispersion Estimation...')
%%% Dispersion estimation %%%
% if isfile(fullfile(savepath,'ProcdData',[fname_save(1:end-7),'_dispCoeff.csv']))
%     dispCoeffs_list = readtable(fullfile(savepath,'CplxProcdData',[fname_save(1:end-7),'_dispCoeff.csv']));
%     dispCoeffs = table2array(dispCoeffs_list);
%     disp('Dispersion coefficient:')
%     disp(dispCoeffs)
% else
%     close all
[dispCoeffs, ref_FFT,ref_FFT_DisComp] = est_dispersion_coeff(fn,LUT,parameters,ref_frame_disp,dispROI,dispMaxOrder);
writematrix(dispCoeffs, fullfile(log_path,[file_id,'_dispCoeff','.csv']));
imwrite(imadjust(mat2gray(20.*log10(abs(ref_FFT(31:1000,:))))),fullfile(log_path,[file_id,'_org.bmp']));
imwrite(imadjust(mat2gray(20.*log10(abs(ref_FFT_DisComp(31:1000,:))))),fullfile(log_path,[file_id,'_disp.bmp']));
disp('Dispersion coefficient:')
disp(dispCoeffs)
% end
% end
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPLEX VOLUME PROCESS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Volume OCT Processing...')
[cplxData_A,cplxData_B] = process_oct_volume(fn,LUT,parameters,dispMaxOrder,dispCoeffs);
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DOPU PROCESS           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('File: %s \n',fname_save);
[OCT_PN, OCT_SN] = est_noise_floor(cplxData_A,cplxData_B,ref_frame_noise);
writematrix([OCT_PN, OCT_SN], fullfile(log_path,[file_id,'_noise_floor.csv']));
ref_noise_FFT = cat(2,cplxData_A(1:1000,:,ref_frame_noise),cplxData_B(1:1000,:,ref_frame_noise));
imwrite(imadjust(mat2gray(20*log10(abs(ref_noise_FFT)))),fullfile(log_path,[file_id,'_noiseFrame.bmp']));

disp('Volume DOPU Processing...')
[avgOCT, DOPU, OCTA] = process_dopu_volume(cplxData_A,cplxData_B, OCT_PN, OCT_SN,numMscans);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save_pdoct_output2(OCT,DOPU,process_path,file_id,cmap_dopu)
save(fullfile(process_path,[file_id,'_avgOCT']), 'avgOCT', '-v7.3');
save(fullfile(process_path,[file_id,'_DOPU']), 'DOPU', '-v7.3');
save(fullfile(process_path,[file_id,'_OCTA']), 'OCTA', '-v7.3');

disp('Saving Complete.');

toc
%%-----------------------------------------------------------------------------------------------------%%
% Preview

% for i = 1:size(avgOCT,3)
%     imagesc(imadjust(mat2gray(abs(cplxData_A(:,:,i))))),colormap('gray')
%     pause(0.1)
% end





