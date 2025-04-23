%% OCT/OCTA Processing Pipeline %% OCT_BatchProcess_VCSEL
% updated : 2022.05.18
clearvars
folder = 'G:\EDUADO\RP005_2022.05.01_OU\OD';
% addpath('C:\Users\KoreanBorg\Desktop\Reseachcore\Project-Melanoma\code');
addpath('C:\Users\coil_\Desktop\Github\Project-PD_OCT\Melanoma\code\OCT_Process');
% folder_LUT = 'C:\Users\KoreanBorg\Desktop\Reseachco       re\Project-Melanoma\code';
folder_LUT = 'C:\Users\coil_\Desktop\PyOCT_YM';
fileIdx = 1;

%%% Preset parameter %%%
dispMaxOrder    = 5;
coeffRange      = 50;
calSigOffSetIdx = 15;

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
fid = fopen(fn);
fseek(fid,2*numPoints*numAscans*(10+round(numBscans/2)-1),-1);
ref_RawData_interlace = fread(fid,[numPoints,numAscans], 'uint16');
ref_RawData = hilbert(cat(2, ref_RawData_interlace(1:2:end,:),ref_RawData_interlace(2:2:end,:)));

% Resampling process %
ref_RawData_rescaled = reSampling_LUT(ref_RawData,LUT);
ref_FFTData_rescaled = fft(ref_RawData_rescaled);
% imagesc(imadjust(mat2gray(20.*log10(abs(ref_FFTData_rescaled(1:end/2,:)))))); colormap(gray)
%% 

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
figure(1),imagesc(imadjust(mat2gray(20.*log10(abs(ref_FFTData_FPNSub(1:600,:)))))); colormap(gray)

%% Dispersion Estimation %%
dispROI        = [150, 400]; % 150, 500 for R005 OS

% Dispersion estimation & compensation %
dispCoeffs_A = setDispCoeff(ref_RawData_HamWin,dispROI,dispMaxOrder,coeffRange);
ref_RawData_DisComp = compDisPhase(ref_RawData_HamWin,dispMaxOrder,dispCoeffs_A);

ref_FFT_Final   = fft(ref_RawData_DisComp);
ref_OCT_Log     = 20.*log10(abs(ref_FFT_Final));
figure(2),imagesc(imadjust(mat2gray(ref_OCT_Log(1:600,:)))); colormap(gray)
%%-----------------------------------------------------------------------------------------------------%%
%% Volume process %%
tic
batchProcessOrder = 1:length(files);
for batchIdx = batchProcessOrder
    clearvars ProcdData
    fn = fnames{batchIdx};
    fid = fopen(fn);
    for FrameNum = 1:numBscans
        fseek(fid,2*numPoints*numAscans*(FrameNum-1),-1);
        rawData_interlace = fread(fid,[numPoints,numAscans], 'uint16');
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
        fprintf('OCT volume %d process: %d\n', batchIdx,FrameNum);
    end
    fclose(fid);
    
    % Save %
    depthROI        = [1, 600];
    cplxData_A = reorder_bscan(ProcdData(depthROI(1):depthROI(2),1:end/2,:),numMscans);
    cplxData_B = reorder_bscan(ProcdData(depthROI(1):depthROI(2),end/2+1:end,:),numMscans);

    [~,fname_save,~] = fileparts(fn);
    exportTiff(flip(20*log10(abs(cplxData_A)+abs(cplxData_B))),[fname_save,'_avgOCT'])
    save(fullfile(cd,[fname_save,'cplxOCT_A']), 'cplxData_A', '-v7.3');
    save(fullfile(cd,[fname_save,'cplxOCT_B']), 'cplxData_B', '-v7.3');
    
    figure(3)
    scan=10+round(numBscans/2)-1;
    fig_cplx=figure(3);
    fig_frame=20.*log10(abs(ProcdData(depthROI(1):depthROI(2),:,scan)));
    imagesc(imadjust(mat2gray(squeeze(fig_frame)))); 
    colormap(gray); title(['Volume: ',fname_save],'Interpreter','none');
    ylabel(['DispROI: ',num2str(dispROI(1)),':',num2str(dispROI(2))]);
    xlabel(['Bscan: ',num2str(scan)]);
    saveas(fig_cplx, fullfile(cd,[fname_save,'fig_cplxOCT.jpg']));
    
end
toc
% save(fullfile(cd,fname_save), 'ProcdData', '-v7.3');


