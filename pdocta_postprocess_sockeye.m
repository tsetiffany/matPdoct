%% OCT/OCTA/DOPU Processing Pipeline %%
% updated : 2021.10.26
%%-----------------------------------------------------------------------------------------------------%%
% % local computer path % %
% datapath = 'I:\04_Sept_MJ_OD_LongRange';
% savepath = datapath;
% fileIdx  = 1;
% addpath(genpath('C:\Users\tiffa\Documents\1. Projects\PDOCT Processing\main_vctrl\Code_PDOCT_OCTA\matPdoct'));

%%% Preset parameter %%%
dispMaxOrder    = 5;
coeffRange      = 50;
bitDepth        = 12;
byteSize        = bitDepth/8;

% cluster path %
SCRIPTLOC = '/scratch/st-mjju-1/tffnytse/Code_PDOCTA/matPdoct';
addpath(genpath(SCRIPTLOC));
%%-----------------------------------------------------------------------------------------------------%%

%%%%% Find filenames of RAW files to be processed %%%%%
cd(datapath);
files   = (dir('*.unp'));
fnames  = {files.name}';
if fileIdx > length(fnames)
    error('Job terminated.\n Index %d exeeds maximum number of files.', fileIdx)
end
fn = fnames{fileIdx};

%%% Load Acquisition Parameters %%%
% Two polarization channels (P/S) are interleaved sample-by-sample within
% each A-scan's raw spectral trace, so numPoints is doubled and numAscans
% is halved relative to the single-channel OCTA convention (matches
% process_oct_volume.m / est_dispersion_coeff.m).
parameters    = getParameters(fn);
numPoints     = parameters(1)*2;
numAscans     = parameters(2)/2;
numBscans     = parameters(3);
numCscans     = parameters(4);
numMscans     = parameters(5);

fprintf('Loaded %s: numPoints=%d, numAscans=%d, numBscans=%d\n', fn, numPoints, numAscans, numBscans);

%%% LUT %%%
fn_ResParam = 'LUT_5050.bin';
fid_ResParam = fopen(fullfile(SCRIPTLOC,fn_ResParam));
rescaleParam = fread(fid_ResParam, 'double');
LUT =  rescaleParam;
fclose(fid_ResParam);

% Create ProcdData output folders %
cd(savepath)
[~,fname_save,~] = fileparts(fn);
pathSplit = regexp(datapath,'\','split');
if length(pathSplit)==1
    pathSplit = regexp(datapath,'/','split');
end
if ~isfolder('ProcdData')
    mkdir ProcdData
end
cd('ProcdData')
if ~isfolder(pathSplit{end-1})
    mkdir(pathSplit{end-1})
end
cd(pathSplit{end-1})
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
log_path     = fullfile(savepath,'ProcdData',pathSplit{end-1},pathSplit{end},fname_save,'log');
file_id      = char(fname_save(end-5:end));

% Output filenames
subOCTA_file     = fullfile(process_path, [file_id '_subOCTA.mat']);
avgOCT_file      = fullfile(process_path, [file_id '_avgOCT.mat']);
nifti_file       = fullfile(process_path, [file_id '_avgOCT']);
dopu_file        = fullfile(process_path, [file_id '_dopu.mat']);

%%-----------------------------------------------------------------------------------------------------%%
%% Reference Frame process %%
ref_frame = 1001;
if bitDepth == 16
    fid = fopen(fn);
    fseek(fid,byteSize*numPoints*numAscans*(ref_frame+1),-1);
    ref_RawData_interlace = fread(fid,[numPoints,numAscans], 'uint16');
elseif bitDepth == 12
    ref_RawData_interlace = unpack_u12u16(fn,numPoints,numAscans,ref_frame);
else
    error('%d-bit data not supported.', bitDepth)
end
ref_RawData = hilbert(cat(2, ref_RawData_interlace(1:2:end,:), ref_RawData_interlace(2:2:end,:)));

% Resampling process %
ref_RawData_rescaled = reSampling_LUT(ref_RawData,LUT);

% FPN removal process (per channel, then recombined) %
ref_RawData_FPNSub_A = ref_RawData_rescaled(:,1:end/2)...
    - (repmat(median(real(ref_RawData_rescaled(:,1:end/2)),2), [1,size(ref_RawData_rescaled(:,1:end/2),2)])...
    +1j.*repmat(median(imag(ref_RawData_rescaled(:,1:end/2)),2), [1,size(ref_RawData_rescaled(:,1:end/2),2)]));
ref_RawData_FPNSub_B = ref_RawData_rescaled(:,end/2+1:end)...
    - (repmat(median(real(ref_RawData_rescaled(:,end/2+1:end)),2), [1,size(ref_RawData_rescaled(:,end/2+1:end),2)])...
    +1j.*repmat(median(imag(ref_RawData_rescaled(:,end/2+1:end)),2), [1,size(ref_RawData_rescaled(:,end/2+1:end),2)]));
ref_RawData_FPNSub = cat(2, ref_RawData_FPNSub_A, ref_RawData_FPNSub_B);

% Windowing process %
ref_RawData_HamWin = ref_RawData_FPNSub...
    .*repmat(hann(size(ref_RawData_FPNSub,1)),[1 size(ref_RawData_FPNSub,2)]);

%% Dispersion Estimation %%
dispROI        = [1, 700];

% Dispersion estimation & compensation (single set, applied to both channels) %
dispCoeffs_A = setDispCoeff(ref_RawData_HamWin,dispROI,dispMaxOrder,coeffRange);
ref_RawData_DisComp = compDisPhase(ref_RawData_HamWin,dispMaxOrder,dispCoeffs_A);

ref_FFT_Final   = fft(ref_RawData_DisComp);
ref_OCT_Log     = 20.*log10(abs(ref_FFT_Final));

clearvars ref_*

disp('Dispersion estimation complete. Starting volume process...');

%% Volume process %%

rawFile = fullfile(datapath, fn);

ProcdData = zeros(1000, numAscans*2, numBscans, 'like', 1i);

for FrameNum = 1:numBscans

    rawData_interlace = unpack_u12u16(rawFile, numPoints, numAscans, FrameNum);
    rawData = hilbert(cat(2, rawData_interlace(1:2:end,:), rawData_interlace(2:2:end,:)));
    rawData_Rescaled = reSampling_LUT(rawData, LUT);

    rawData_FPNSub_A = rawData_Rescaled(:,1:end/2) - ...
        (repmat(median(real(rawData_Rescaled(:,1:end/2)),2), [1 size(rawData_Rescaled(:,1:end/2),2)]) + ...
         1j .* repmat(median(imag(rawData_Rescaled(:,1:end/2)),2), [1 size(rawData_Rescaled(:,1:end/2),2)]));
    rawData_FPNSub_B = rawData_Rescaled(:,end/2+1:end) - ...
        (repmat(median(real(rawData_Rescaled(:,end/2+1:end)),2), [1 size(rawData_Rescaled(:,end/2+1:end),2)]) + ...
         1j .* repmat(median(imag(rawData_Rescaled(:,end/2+1:end)),2), [1 size(rawData_Rescaled(:,end/2+1:end),2)]));
    rawData_FPNSub = cat(2, rawData_FPNSub_A, rawData_FPNSub_B);

    rawData_HamWin = rawData_FPNSub .* repmat(hann(size(rawData_FPNSub,1)), [1 size(rawData_FPNSub,2)]);

    rawData_DisComp = compDisPhase(rawData_HamWin, dispMaxOrder, dispCoeffs_A);
    fftData_DispComp = fft(rawData_DisComp);

    ProcdData(:,:,FrameNum) = fftData_DispComp(1:1000,:);

end

clearvars fftData_DispComp rawData*

disp('Volume process complete. Splitting polarization channels...');

%% Split polarization channels %%
% ProcdData columns 1:end/2 are channel A (P), end/2+1:end are channel B (S)

cplxData_A = ProcdData(:,1:end/2,:);
cplxData_B = ProcdData(:,end/2+1:end,:);
clear ProcdData

disp('Channels split. Starting interleave...');

%% Interleave (per channel) %%
% The forward/backward 4BM-scan interleave + lateral crop (51:550)
% Step-bidirectional (0101 0101 2323 2323...)

cplxData_A = interleaveChannel(cplxData_A);
cplxData_B = interleaveChannel(cplxData_B);

disp('Interleave complete.');

%% Check Parameters
numPoints  = size(cplxData_A,1);
numAscans  = size(cplxData_A,2);
numBscans  = size(cplxData_A,3);
numBMscans = 4; % Number of repeated B-scans

%% DOPU Processing %%
ref_frame_noise = 2;

[OCT_PN, OCT_SN] = est_noise_floor(cplxData_A,cplxData_B,ref_frame_noise);
writematrix([OCT_PN, OCT_SN], fullfile(log_path,[file_id,'_noise_floor.csv']));
ref_noise_FFT = cat(2,cplxData_A(:,:,ref_frame_noise),cplxData_B(:,:,ref_frame_noise));
imwrite(imadjust(mat2gray(20*log10(abs(ref_noise_FFT)))),fullfile(log_path,[file_id,'_noiseFrame.bmp']));
clear ref_noise_FFT

disp('Starting DOPU processing...');
[DOPU] = process_dopu_volume(cplxData_A, cplxData_B, OCT_PN, OCT_SN, numBMscans);

%% OCTA Processing %%
disp('Starting OCTA (channel A)...');

usfac = 20; % Up-sampling Factor

[avgOCT_A, subOCTA_A] = process_octa_volume(cplxData_A, numBMscans, usfac);
clear cplxData_A

disp('Starting OCTA (channel B)...');

[avgOCT_B, subOCTA_B] = process_octa_volume(cplxData_B, numBMscans, usfac);
clear cplxData_B

disp('Averaging channels...');

avgOCT  = (avgOCT_A + avgOCT_B) / 2;
subOCTA = (subOCTA_A + subOCTA_B) / 2;

clear avgOCT_A avgOCT_B subOCTA_A subOCTA_B

disp('OCTA calculation complete. Saving outputs...');

%%
save(subOCTA_file, "subOCTA", '-v7.3');
save(avgOCT_file, "avgOCT", '-v7.3');
save(dopu_file, "DOPU", '-v7.3');
clear DOPU

imadjusted_vol = zeros(size(avgOCT));
for i = 1:size(avgOCT,3)
    imadjusted_vol(:,:,i) = imadjust(mat2gray(avgOCT(:,:,i)));
end

niftiwrite(imadjusted_vol, nifti_file, 'Compressed', true);

disp('Saving complete.');

%% Local functions %%

function out = interleaveChannel(vol)
% Reproduces the forward/backward 4BM-scan interleave and 51:550 lateral
% crop (same logic as the single-channel OCTA pipeline), operating on
% one channel's own A-scan width.

    fwdScan = vol(:,:,1:2:end);
    bwdScan = fliplr(vol(:,:,2:2:end));

    fwdScan = fwdScan(:,51:550,:);
    bwdScan = bwdScan(:,51:550,:);

    blockSize = 4;
    numFrames = size(fwdScan,3);

    out = zeros(size(fwdScan,1), size(fwdScan,2), size(fwdScan,3)+size(bwdScan,3), 'like', fwdScan);

    outIdx = 1;
    for k = 1:blockSize:numFrames

        % Forward block
        out(:,:,outIdx:outIdx+blockSize-1) = fwdScan(:,:,k:k+blockSize-1);
        outIdx = outIdx + blockSize;

        % Backward block
        out(:,:,outIdx:outIdx+blockSize-1) = bwdScan(:,:,k:k+blockSize-1);
        outIdx = outIdx + blockSize;

    end
end
