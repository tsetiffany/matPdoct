clearvars
% folder = 'D:\Data\2021.09.21\calibration';
folder = 'D:\Data\2021.10.05\test';
addpath('C:\Users\KoreanBorg\Desktop\Reseachcore\Matlab-FDML_Postprocessing');
addpath('C:\Users\KoreanBorg\Desktop\Reseachcore\Matlab-VCSEL_Postprocessing');

cd(folder);
files   = (dir('*.unp'));
fnames  = {files.name}';

%%% LUT %%%
fn_ResParam  = 'LUTSS.bin';
fid_ResParam = fopen(fn_ResParam);
LUT = fread(fid_ResParam, 'double');
fclose(fid_ResParam);

%%% Load Data %%%
ref_Frame = 17; % round(numBscans/2);
i=1;
% for i = 1:length(fnames)
    fn = fnames{i};
    parameters    = getParameters(fn);
    numPoints     = parameters(1);
    numAscans     = parameters(2);
    numBscans     = parameters(3);
    numCscans     = parameters(4);
    numMscans     = parameters(5);

    fid = fopen(fn);
    fseek(fid,2*numPoints*numAscans*(ref_Frame-1),-1);
    rawData  = hilbert(fread(fid,[numPoints,numAscans], 'uint16'));
%     rawData_rescaled(:,:,i) = reSampling_LUT(rawData, LUT);
    fclose(fid);
% end
% dispMaxOrder    = 5;
% dispCoeffs_A = [3.3334,0.1352,0.0292,0.0832];

%% Spectrum Test %%
% ref_data = 1;

ref_rawData = rawData;

figure,plot(real(ref_rawData(2:2:end,252))); 
% xlim([-50 1280])
ylim([0 55000])

for i=1:4
    subplot(2,2,i),plot(real(ref_rawData(:,i+199))),title(sprintf('frame=%d',i+199))
end

ref_fftData   = fft(ref_rawData);

plot(20.*log10(abs(mean(ref_fftData(1:500,:),2))));

calSigOffSetIdx = 15;
diffPhasesA = estPhaseShift(ref_fftData(:,1), ref_fftData, calSigOffSetIdx);
figure(1),
subplot(1,2,1),plot(diffPhasesA)
subplot(1,2,2),histogram(diffPhasesA)

Win_cal    = [59, 116];
Win_sig    = [251, 305];
winFunc_cal = zeros(size(ref_fftData));
winFunc_cal(Win_cal(1):Win_cal(2),:) = 1;
winFunc_sig = zeros(size(ref_fftData));
winFunc_sig(Win_sig(1):Win_sig(2),:) = 1;

cal_fftData_extr = ref_fftData.*winFunc_cal;
sig_fftData_extr = ref_fftData.*winFunc_sig;

rawData_extr_cal = ifft(cal_fftData_extr);
rawData_extr_sig = ifft(sig_fftData_extr);
plot(real(mean(rawData_extr_cal,2))), hold on
plot(real(mean(rawData_extr_sig,2))), hold off

diffPhasesA = estPhaseDiff(rawData_extr_cal(:,1), rawData_extr_cal, calSigOffSetIdx);
diffPhasesB = estPhaseDiff(rawData_extr_sig(:,1), rawData_extr_sig, calSigOffSetIdx);
figure(1),
subplot(1,2,1),plot(diffPhasesA),title('calibration'),ylim([-3.5 3.5])
subplot(1,2,2),plot(diffPhasesB),title('signal'),ylim([-3.5 3.5])

spectrum_1 = abs(rawData_extr);

spectrum_1 = spectrum_1./(max(spectrum_1(:)));

plot(spectrum_1); 
set(gca,'FontSize',17);

phase_1 = unwrap(angle(rawData_extr));
phase_1 = phase_1/phase_1(size(phase_1,1)/2);
phase_1_norm = phase_1 - min(phase_1(:));
phase_1_norm = phase_1_norm ./ max(phase_1_norm(:));

plot(phase_1_norm); 

%% apply LUT %%

rawData_avg   = squeeze(mean(rawData,2));
fftData_avg   = fft(rawData_avg);

plot(20.*log10(abs(fftData_avg(1:end/2,:))))

rawDataArray = rawData_avg;
rawDataArray_rescaled = reSampling_LUT(rawDataArray, LUT);

rawDataArray_han = rawDataArray...
    .*repmat(hann(size(rawDataArray,1)),[1 size(rawDataArray,2)]);

rawDataArray_rescaled_han = rawDataArray_rescaled...
    .*repmat(hann(size(rawDataArray_rescaled,1)),[1 size(rawDataArray_rescaled,2)]);

rawDataArray_dispComp_han = compDisPhase(rawDataArray_rescaled_han,dispMaxOrder,dispCoeffs_A);

% fftDataArray = fft(rawDataArray_han);
% fftDataArray_rescaled = fft(rawDataArray_rescaled_han);
fftDataArray = fft([rawDataArray_han; zeros(size(rawDataArray,1)*9,size(rawDataArray,2))]);
fftDataArray_rescaled = fft([rawDataArray_rescaled_han; zeros(size(rawDataArray,1)*9,size(rawDataArray,2))]);
fftDataArray_dispComp = fft([rawDataArray_dispComp_han; zeros(size(rawDataArray,1)*9,size(rawDataArray,2))]);

subplot(1,3,1)
plot(20.*log10(abs(fftDataArray(1:end/2,:)))-80)
ylim([-10 80])
title('Before k-lienar')
subplot(1,3,2)
plot(20.*log10(abs(fftDataArray_rescaled(1:end/2,:)))-60)
ylim([-10 80])
title('After k-linear')
subplot(1,3,3)
plot(20.*log10(abs(fftDataArray_dispComp(1:end/2,:)))-60)
ylim([-10 80])
title('After k-linear & DC')

%% Estimate dispersion compensation %% 

rawDataArray_rescaled = rawData_rescaled;

skip = 16;
rawDataArray_cat  = cat(2,rawDataArray_rescaled(:,1:skip:end,1),rawDataArray_rescaled(:,1:skip:end,2),...
    rawDataArray_rescaled(:,1:skip:end,3),rawDataArray_rescaled(:,1:skip:end,4),...
    rawDataArray_rescaled(:,1:skip:end,5),rawDataArray_rescaled(:,1:skip:end,6),...
    rawDataArray_rescaled(:,1:skip:end,7),rawDataArray_rescaled(:,1:skip:end,8),...
    rawDataArray_rescaled(:,1:skip:end,9),rawDataArray_rescaled(:,1:skip:end,10),...
    rawDataArray_rescaled(:,1:skip:end,11),rawDataArray_rescaled(:,1:skip:end,12),...
    rawDataArray_rescaled(:,1:skip:end,13),rawDataArray_rescaled(:,1:skip:end,14));
rawDataArray_han = rawDataArray_cat...
    .*repmat(hann(size(rawDataArray_cat,1)),[1 size(rawDataArray_cat,2)]);
fftDataArray_han = fft(rawDataArray_han);
figure(1), imagesc(imadjust(mat2gray(abs(fftDataArray_han(1:end/2,:))))); colormap(gray)

% Dispersion estimation & compensation %
dispMaxOrder    = 5;
coeffRange      = 50;
depthROI        = [80, 600];

% dispCoeffs_A = setDispCoeff(rawDataArray_han,depthROI,dispMaxOrder,coeffRange);
ref_FFTData_DisComp = fft(compDisPhase(rawDataArray_han,dispMaxOrder,dispCoeffs_A));
figure(2), imagesc(imadjust(mat2gray(abs(ref_FFTData_DisComp(1:end/2,:))))); colormap(gray)

% % Save
% fileID = fopen('dispCoeffs.bin','w');
% fwrite(fileID,dispCoeffs_A,'double');
% fclose(fileID);



