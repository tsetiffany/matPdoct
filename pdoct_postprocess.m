%% OCT/OCTA Processing Pipeline %%
% updated : 2025.03.17

%%
%%-----------------------------------------------------------------------------------------------------%%
% local computer path % %
tic
datapath = 'I:\Ansel_OCTA_500\test';
code_loc = 'C:\Users\tiffa\Documents\1. Projects\PDOCT Processing\main_vctrl\Code_workingversion\matPdoct';
addpath(genpath(code_loc));
fileIdx = 1;
savepath = datapath;
dispROI = [30, 500];

% cluster path %
% SCRIPTLOC = '/scratch/st-mjju-1/tffnytse/Code/matPdoct';
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
% fn_ResParam = 'LUTSS.bin';
fn_ResParam = 'LUT_5050.bin';
fid_ResParam = fopen(fn_ResParam);
rescaleParam = fread(fid_ResParam, 'double');
LUT =  rescaleParam;
fclose(fid_ResParam);

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

[dispCoeffs, ref_FFT,ref_FFT_DisComp] = est_dispersion_coeff(fn,LUT,parameters,ref_frame_disp,dispROI,dispMaxOrder);
writematrix(dispCoeffs, fullfile(log_path,[file_id,'_dispCoeff','.csv']));
imwrite(imadjust(mat2gray(20.*log10(abs(ref_FFT(31:1000,:))))),fullfile(log_path,[file_id,'_org.bmp']));
imwrite(imadjust(mat2gray(20.*log10(abs(ref_FFT_DisComp(31:1000,:))))),fullfile(log_path,[file_id,'_disp.bmp']));

disp('Dispersion coefficient:')

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPLEX VOLUME PROCESS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('Volume OCT Processing...')
[cplxData_A,cplxData_B] = process_oct_volume(fn,LUT,parameters,dispMaxOrder,dispCoeffs);

cplxData_A = reorder_bscan(cplxData_A,numMscans);
cplxData_B = reorder_bscan(cplxData_B,numMscans);
save(fullfile(process_path,[file_id,'_cplxData_A']),'cplxData_A','-v7.3')
save(fullfile(process_path,[file_id,'_cplxData_B']),'cplxData_B','-v7.3')

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OCTA PROCESS           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[avgOCT, avgOCT_tcorr, OCTA, OCTA_tcorr] = process_octa_volume(cplxData_A,cplxData_B,numMscans,process_path);


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(fullfile(process_path,[file_id,'_avgOCT']), 'avgOCT', '-v7.3');
save(fullfile(process_path,[file_id,'_avgOCT_tcorr']), 'avgOCT_tcorr', '-v7.3');
save(fullfile(process_path,[file_id,'_OCTA']), 'OCTA', '-v7.3');
save(fullfile(process_path,[file_id,'_OCTA_tcorr']), 'OCTA_tcorr', '-v7.3');

disp('Saving Complete.');

toc




