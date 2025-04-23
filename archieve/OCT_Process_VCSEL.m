%% OCT/OCTA Processing Pipeline %%
% updated : 2021.10.26
tic
clearvars, clc, close all
folder = 'G:\VCSEL_PDOCT\2022.12.14\MIS';
addpath(genpath('C:\Users\coil_\OneDrive\Documents\Github\yusi-miao\Software-MatOctPipeline\PDOCT_400VCSEL'));
folder_LUT = 'C:\Users\coil_\OneDrive\Documents\Github\yusi-miao\Software-MatOctPipeline\PDOCT_400VCSEL';
fileIdx = 2;

%%% Preset parameter %%%
dispMaxOrder    = 5;
coeffRange      = 50;
calSigOffSetIdx = 15;
bitDepth        = 12;
byteSize        = bitDepth/8;

%%%%% Find filenames of RAW files to be processed %%%%%
cd(folder);
files   = (dir('*.unp'));
fnames  = {files.name}';
fn = fnames{fileIdx};

%%% Load Acquisition Parameters %%%
parameters    = getParameters(fn);
numPoints     = parameters(1)*2;
numAscans     = parameters(2)/2;
numBscans     = parameters(3);
numCscans     = parameters(4);
numMscans     = parameters(5);

%%% LUT %%%
fn_ResParam = 'LUTSS.bin';
fid_ResParam = fopen(fullfile(folder_LUT,fn_ResParam));
rescaleParam = fread(fid_ResParam, 'double');
LUT =  rescaleParam;
fclose(fid_ResParam);

%%-----------------------------------------------------------------------------------------------------%%
%% Reference Frame process %%
% ref_frame = round(numBscans/2);
ref_frame = 511;
if bitDepth == 16
    fid = fopen(fn);
    fseek(fid,byteSize*numPoints*numAscans*(ref_frame+1),-1);
    ref_RawData_interlace = fread(fid,[numPoints,numAscans], 'uint16');
elseif bitDepth == 12
    ref_RawData_interlace = unpack_u12u16(fn,numPoints,numAscans,ref_frame);
else
    error('%d-bit data not supported.', bitDepth)
end
ref_RawData = hilbert(cat(2, ref_RawData_interlace(1:2:end,:),ref_RawData_interlace(2:2:end,:)));
figure(1),plot(real(ref_RawData(:,1))), hold on, plot(real(ref_RawData(:,end/2+1))),hold off,

% Resampling process %
ref_RawData_rescaled = reSampling_LUT(ref_RawData,LUT);
ref_FFTData_rescaled = fft(ref_RawData_rescaled);
% imagesc(imadjust(mat2gray(20.*log10(abs(ref_FFTData_rescaled(1:end/2,:)))))); colormap(gray)

% FPN removal process %
ref_RawData_FPNSub_A = ref_RawData_rescaled(:,1:end/2)...
    - (repmat(median(real(ref_RawData_rescaled(:,1:end/2)),2), [1,size(ref_RawData_rescaled(:,1:end/2),2)])...
    +1j.*repmat(median(imag(ref_RawData_rescaled(:,1:end/2)),2), [1,size(ref_RawData_rescaled(:,1:end/2),2)]));
ref_RawData_FPNSub_B = ref_RawData_rescaled(:,end/2+1:end)...
    - (repmat(median(real(ref_RawData_rescaled(:,end/2+1:end)),2), [1,size(ref_RawData_rescaled(:,end/2+1:end),2)])...
    +1j.*repmat(median(imag(ref_RawData_rescaled(:,end/2+1:end)),2), [1,size(ref_RawData_rescaled(:,end/2+1:end),2)]));
ref_RawData_FPNSub = cat(2,ref_RawData_FPNSub_A,ref_RawData_FPNSub_B);
ref_FFTData_FPNSub = fft(ref_RawData_FPNSub);

% Windowing process %
ref_RawData_HamWin = ref_RawData_FPNSub...
    .*repmat(hann(size(ref_RawData_FPNSub,1)),[1 size(ref_RawData_FPNSub,2)]);
ref_FFTData_HamWin = fft(ref_RawData_HamWin);
figure(1),imagesc(imadjust(mat2gray(20.*log10(abs(ref_FFTData_FPNSub(51:750,:)))))); colormap(gray)

%% Dispersion Estimation %%
dispROI        = [100, 500];

% Dispersion estimation & compensation %
dispCoeffs_A = setDispCoeff(ref_RawData_HamWin,dispROI,dispMaxOrder,coeffRange);
ref_RawData_DisComp = compDisPhase(ref_RawData_HamWin,dispMaxOrder,dispCoeffs_A);

ref_FFT_Final   = fft(ref_RawData_DisComp);
ref_OCT_Log     = 20.*log10(abs(ref_FFT_Final));
figure(2),imagesc(imadjust(mat2gray(ref_OCT_Log(51:750,:)))); colormap(gray)
%%-----------------------------------------------------------------------------------------------------%%
toc
%% Volume process %%
tic
clearvars ProcdData
for FrameNum = 1:numBscans
    fseek(fid,2*numPoints*numAscans*(FrameNum-1),-1);
    rawData_interlace = fread(fid,[numPoints,numAscans], 'ubit12');
    rawData = hilbert(cat(2, rawData_interlace(1:2:end,:),rawData_interlace(2:2:end,:)));
    
    % Rescale process %
    rawData_Rescaled = reSampling_LUT(rawData,LUT);
    fftData_Rescaled = fft(rawData_Rescaled);
    
    % FPN remove %
    rawData_FPNSub_A = rawData_Rescaled(:,1:end/2)...
        - (repmat(median(real(rawData_Rescaled(:,1:end/2)),2), [1,size(rawData_Rescaled(:,1:end/2),2)])...
        +1j.*repmat(median(imag(rawData_Rescaled(:,1:end/2)),2), [1,size(rawData_Rescaled(:,1:end/2),2)]));
    rawData_FPNSub_B = rawData_Rescaled(:,end/2+1:end)...
        - (repmat(median(real(rawData_Rescaled(:,end/2+1:end)),2), [1,size(rawData_Rescaled(:,end/2+1:end),2)])...
        +1j.*repmat(median(imag(rawData_Rescaled(:,end/2+1:end)),2), [1,size(rawData_Rescaled(:,end/2+1:end),2)]));
    rawData_FPNSub = cat(2,rawData_FPNSub_A,rawData_FPNSub_B);
    
    % Windowing process %
    rawData_HamWin = rawData_FPNSub...
        .*repmat(hann(size(rawData_FPNSub,1)),[1 size(rawData_FPNSub,2)]);

    % Dispersion estimation & compensation %
    rawData_DisComp = compDisPhase(rawData_HamWin,dispMaxOrder,dispCoeffs_A);
    fftData_DispComp = fft(rawData_DisComp);

    ProcdData(:,:,FrameNum) = fftData_DispComp(1:1100,:);
    fprintf('OCT volume process: %d\n', FrameNum);
end
toc
% save(fullfile(cd,fname_save), 'ProcdData', '-v7.3');

% Save %
depthROI        = [1, 1000];
cplxData_A = reorder_bscan(ProcdData(depthROI(1):depthROI(2),1:end/2,:),numMscans);
cplxData_B = reorder_bscan(ProcdData(depthROI(1):depthROI(2),end/2+1:end,:),numMscans);

[~,fname_save,~] = fileparts(fn);
exportTiff(flip(20*log10(abs(cplxData_A)+abs(cplxData_B))),[fname_save,'_avgOCT'])
% exportTiff(flip(20*log10(abs(cplxData_B))),[fname_save,'OCT_B'])
save(fullfile(cd,[fname_save,'A']), 'cplxData_A', '-v7.3');
save(fullfile(cd,[fname_save,'B']), 'cplxData_B', '-v7.3');
