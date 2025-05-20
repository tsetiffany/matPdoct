folder = 'C:\Users\mjju\Desktop\SAO-OCTA';
addpath('C:\Users\mjju\Dropbox\ProgramScripts\MatlabScripts\OCTViewer_Project');

rescaleFolder = 'C:\Users\mjju\Desktop\SAO-OCTA';
fn_ResParam = 'LUTSS.dat';
fid_ResParam = fopen(fullfile(rescaleFolder,fn_ResParam));
rescaleParam = fread(fid_ResParam, 'double')+1.00;
LUT =  rescaleParam;
fclose(fid_ResParam);

dispMaxOrder    = 5;
coeffRange      = 20;

%%%%% Find filenames of RAW files to be processed %%%%%
cd(folder);
files   = (dir('*.unp'));
fnames  = {files.name}';
fn = fnames{1};

%%% Load Acquisition Parameters %%%
parameters    = getParameters(fn);
numPoints     = parameters(1);
numAscans     = parameters(2);
numBscans     = parameters(3);
numCscans     = parameters(4);
numMscans     = parameters(5);
calSignal     = [1 25]; %[idx_start, idx_end]

fid = fopen(fn);


winFunc_Img     = ones([numPoints, 1]);
% winFunc_Img(1:calSignal(2)+10) = 0;
winFunc_Cal     = zeros([numPoints, 1]);
winFunc_Cal(calSignal(1):calSignal(2)) = 1;
winFunc_Cal2     = zeros([numPoints, 1]);
winFunc_Cal2(15:calSignal(2)) = 1;

ref_Frame        = 100; % select reference frame for estimating dispersion coefficients.
fseek(fid,2*numPoints*numAscans*(ref_Frame-1),-1);
ref_RawData      = fread(fid,[numPoints,numAscans], 'uint16');
ref_FFTData      = fft(hilbert(ref_RawData));

ref_FFTData_Img  = ref_FFTData;
ref_FFTData_Cal  = ref_FFTData(:,100).*winFunc_Cal;

ref_RawData_Img  = ifft(ref_FFTData_Img);
ref_RawData_Cal  = ifft(ref_FFTData_Cal);
ref_Img      = reSampling(fft(ref_RawData_Img), LUT);
img_FPNSub   = ref_Img - repmat(mean(ref_Img,2), [1,size(ref_Img,2)]);
ref_Img_HamWin   = img_FPNSub.*repmat(hamming(size(ref_Img,1)),[1 size(ref_Img,2)]);
dispCoeffs   = setDispCoeff(ref_Img_HamWin, dispMaxOrder, coeffRange);

ref_DispComp     = compDisPhase(ref_Img_HamWin,dispMaxOrder,dispCoeffs);
ref_FFT_Img = fft(ref_DispComp);
ref_OCT = 20.*log(abs(ref_FFT_Img(1:400,:)));
figure
imagesc(mat2gray(ref_OCT)); colormap(gray)

plot(real(ref_Img_HamWin))

for FrameNum = 1:numBscans    
    fseek(fid,2*numPoints*numAscans*(FrameNum-1),-1);
    
    %%% Load raw spectrum %%%
    rawData       = fread(fid,[numPoints,numAscans], 'uint16');
    fftData       = fft(hilbert(rawData));
    
%     fftData_Img   = fftData.*repmat(winFunc_Img,[1 size(fftData,2)]);
    fftData_Img   = fftData;
    fftData_Cal   = fftData.*repmat(winFunc_Cal,[1 size(fftData,2)]);
    fftData_Cal2  = fftData.*repmat(winFunc_Cal2,[1 size(fftData,2)]);
    
    %%% Estimating spectral-shift amount in sub-pixel %%%
    phaseDiff    = estPhaseDiff2(ref_FFTData_Cal, fftData_Cal2);
    phaseSlop    = getDiffPhase(fftData_Cal2, LUT, phaseDiff);
    
    %%% Spectral shift compensation during resampling process %%%
    img           = reSampling(fftData_Cal, LUT);
    img_phaseComp = ifft(fft(img).*exp(-1j.*phaseSlop));

    %%% DC & Fixed pattern noise remove process through medial filtering %%%
    img_FPNSub    = img_phaseComp...
        - (repmat(median(real(img_phaseComp),2), [1,size(img_phaseComp,2)])...
        +1j.*repmat(median(imag(img_phaseComp),2), [1,size(img_phaseComp,2)]));

    %%% Dispersion compensation with the estimated dispersion coefficients %%%
    img_DispComp  = compDisPhase(img_FPNSub,dispMaxOrder,dispCoeffs);
    
    %%% Windowing process %%%
    img_HamWin    = img_DispComp.*repmat(hann(size(img_DispComp,1)),[1 size(img_DispComp,2)]);
    
    %%% FFT %%%
    img_FFT                 = fft(img_HamWin);

    ProcdData(:,:,FrameNum) = img_FFT(1:400,:,:);
end




