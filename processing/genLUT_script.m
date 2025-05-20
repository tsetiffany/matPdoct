folder = 'C:\Rodent\RawData\2021.02.09\570';
% folder = 'C:\Rodent\RawData\2021.02.16\CalData_570';
addpath('C:\Users\KoreanBorg\Desktop\Reseachcore\Matlab-FDML_Postprocessing');

cd(folder);
files   = (dir('*.unp'));
fnames  = {files.name}';

for i = 1:length(fnames)
    fn = fnames{i};
    parameters    = getParameters(fn);
    numPoints     = parameters(1);
    numAscans     = parameters(2);
    numBscans     = parameters(3);
    numCscans     = parameters(4);
    numMscans     = parameters(5);

    fid = fopen(fn);
    fseek(fid,2*numPoints*numAscans*((numBscans/2)-1),-1);
    rawData(:,:,i)  = hilbert(fread(fid,[numPoints,numAscans], 'uint16'));
end

%% Spectrum Test %%
ref_data = 6;

rawData_1   = squeeze(mean(rawData(:,1:4:end,ref_data),2));
rawData_2   = squeeze(mean(rawData(:,2:4:end,ref_data),2));
rawData_3   = squeeze(mean(rawData(:,3:4:end,ref_data),2));
rawData_4   = squeeze(mean(rawData(:,4:4:end,ref_data),2));

plot(real(rawData_1)); hold on
% plot(real(rawData_2));
% plot(real(rawData_3));
plot(real(rawData_4)); hold off

% offSet = [1, 1152];
offSet = [25, 1134];

fftData_1   = fft(rawData_1(offSet(1):offSet(2),:));
fftData_2   = fft(rawData_2(offSet(1):offSet(2),:));
fftData_3   = fft(rawData_3(offSet(1):offSet(2),:));
fftData_4   = fft(rawData_4(offSet(1):offSet(2),:));

plot(20.*log10(abs(fftData_1(1:500,:)))); hold on
plot(20.*log10(abs(fftData_2(1:500,:))));
plot(20.*log10(abs(fftData_3(1:500,:)))); 
plot(20.*log10(abs(fftData_4(1:500,:)))); hold off

Win    = [152, 243];
% Win    = [175, 278];
% Win    = [199, 276];
winFunc = zeros(size(fftData_1));
winFunc(Win(1):Win(2),:) = 1;

fftData_extr_1 = fftData_1.*winFunc;
fftData_extr_2 = fftData_2.*winFunc;
fftData_extr_3 = fftData_3.*winFunc;
fftData_extr_4 = fftData_4.*winFunc;

rawData_extr_1 = ifft(fftData_extr_1);
rawData_extr_2 = ifft(fftData_extr_2);
rawData_extr_3 = ifft(fftData_extr_3);
rawData_extr_4 = ifft(fftData_extr_4);

spectrum_1 = abs(rawData_extr_1);
spectrum_2 = abs(rawData_extr_2);
spectrum_3 = abs(rawData_extr_3);
spectrum_4 = abs(rawData_extr_4);

spectrum_1 = spectrum_1./(max(spectrum_1(:)));
spectrum_2 = spectrum_2./(max(spectrum_2(:)));
spectrum_3 = spectrum_3./(max(spectrum_3(:)));
spectrum_4 = spectrum_4./(max(spectrum_4(:)));

plot(spectrum_1); hold on
plot(spectrum_2); grid on
plot(spectrum_3);
plot(spectrum_4); hold off
set(gca,'FontSize',17);

phase_1 = unwrap(angle(rawData_extr_1));
phase_2 = unwrap(angle(rawData_extr_2));
phase_3 = unwrap(angle(rawData_extr_3));
phase_4 = unwrap(angle(rawData_extr_4));

phase_1 = phase_1/phase_1(size(phase_1,1)/2);
phase_1_norm = phase_1 - min(phase_1(:));
phase_1_norm = phase_1_norm ./ max(phase_1_norm(:));

phase_2 = phase_2/phase_2(size(phase_2,1)/2);
phase_2_norm = phase_2 - min(phase_2(:));
phase_2_norm = phase_2_norm ./ max(phase_2_norm(:));

phase_3 = phase_3/phase_3(size(phase_3,1)/2);
phase_3_norm = phase_3 - min(phase_3(:));
phase_3_norm = phase_3_norm ./ max(phase_3_norm(:));

phase_4 = phase_4/phase_4(size(phase_4,1)/2);
phase_4_norm = phase_4 - min(phase_4(:));
phase_4_norm = phase_4_norm ./ max(phase_4_norm(:));

plot(phase_1_norm); hold on
plot(phase_2_norm); grid on
plot(phase_3_norm);
plot(phase_4_norm); hold off

phaseDiff12 = phase_1_norm - phase_2_norm;
phaseDiff13 = phase_1_norm - phase_3_norm;
phaseDiff14 = phase_1_norm - phase_4_norm;

% phaseDiff12 = phase_1 - phase_2;
% phaseDiff13 = phase_1 - phase_3;
% phaseDiff14 = phase_1 - phase_4;

plot(phaseDiff12); hold on;
plot(phaseDiff13);
plot(phaseDiff14); hold off;
set(gca,'FontSize',17);
legend({'Phase1 - Phase2','Phase1 - Phase3','Phase1 - Phase4'},'Location','north');

%% LUT generate %%

%First Buffer%
rawData_avg   = squeeze(mean(rawData(offSet(1):offSet(2),4:4:end,:),2));
fftData_avg   = fft(rawData_avg);

plot(20.*log10(abs(fftData_avg(1:end/2,:))))


rawData_cal1 = rawData_avg(:,6);
rawData_cal2 = rawData_avg(:,4);

fftData_cal1 = fft(rawData_cal1);
fftData_cal2 = fft(rawData_cal2);

calWin1    = [170, 278];
% calWin1      = [150, 248];
% calWin1    = [199, 276];
winFunc1 = zeros(size(fftData_cal1));
winFunc1(calWin1(1):calWin1(2),:) = 1;
% win1 = hann(calWin1(2)-calWin1(1));
% win1 = conv(winFunc1(:,1), win1, 'same');
% winFunc1 = repmat(win1, [1,size(winFunc1,2)]);
calFFTData1 = fftData_cal1.*winFunc1;
plot(20.*log10(abs(fftData_cal1(1:end/2,:)))); hold on; plot(20.*log10(abs(calFFTData1(1:end/2,:)))); hold off

calWin2    = [88, 155];
% calWin2    = [64, 123];
% calWin2    = [107, 153];
winFunc2 = zeros(size(fftData_cal2));
winFunc2(calWin2(1):calWin2(2),:) = 1;
% win2 = hann(calWin2(2)-calWin2(1));
% win2 = conv(winFunc2(:,1), win2, 'same');
% winFunc2 = repmat(win2, [1,size(winFunc2,2)]);
calFFTData2 = fftData_cal2.*winFunc2;
plot(20.*log10(abs(fftData_cal2(1:end/2,:)))); hold on; plot(20.*log10(abs(calFFTData2(1:end/2,:)))); hold off


calRawData1 = ifft(calFFTData1);
calRawData2 = ifft(calFFTData2);

calRawData  = calRawData1.*conj(calRawData2);

subplot(1,2,1)
plot(abs(calRawData));
subplot(1,2,2)
plot(unwrap(angle(calRawData)));

%m = (1:size(calRawData,1))';
m = (1:size(calRawData,1))';

q_raw = unwrap(angle(calRawData));

q = q_raw/q_raw(size(q_raw,1)/2);
q_norm = q - min(q(:));
q_norm = q_norm ./ max(q_norm(:));

% p = fit(q_norm, m,'cubicinterp');
p  = polyfit(q_norm, m, 12);

q_lin = linspace(0, 1,size(calRawData,1))';

% newIdx = p(q_lin);
newIdx = polyval(p, q_lin);

LUT = newIdx - min(newIdx(:));
LUT = LUT ./ max(LUT(:));

LUT_acuqi = (LUT.*(size(calRawData,1)-1));
LUT_procd = (LUT.*(size(calRawData,1)));

% fileID = fopen('LUTSS.bin','w');
% fwrite(fileID,LUT_acuqi,'double');
% fclose(fileID);
% 
% fileID = fopen('LUTSS_Procd.bin','w');
% fwrite(fileID,LUT_procd,'double');
% fclose(fileID);


rawDataArray = rawData_avg;
rawDataArray_rescaled = reSampling_LUT(rawDataArray, LUT_procd);

rawDataArray_han = rawDataArray...
    .*repmat(hann(size(rawDataArray,1)),[1 size(rawDataArray,2)]);

rawDataArray_rescaled_han = rawDataArray_rescaled...
    .*repmat(hann(size(rawDataArray_rescaled,1)),[1 size(rawDataArray_rescaled,2)]);

% fftDataArray = fft(rawDataArray_han);
% fftDataArray_rescaled = fft(rawDataArray_rescaled_han);
fftDataArray = fft([rawDataArray; zeros(size(rawDataArray,1)*9,size(rawDataArray,2))]);
fftDataArray_rescaled = fft([rawDataArray_rescaled; zeros(size(rawDataArray,1)*9,size(rawDataArray,2))]);

subplot(1,2,1)
plot(20.*log10(abs(fftDataArray(1:end/2,:)))-80)
% ylim([60 160])
subplot(1,2,2)
plot(20.*log10(abs(fftDataArray_rescaled(1:end/2,:)))-60)
% ylim([40 120])

%% Spectral testing after resampling

rawData_1_resample = reSampling_LUT(rawData_1, LUT_procd);
% rawData_2_resample = reSampling_LUT(rawData_2, LUT_procd);
% rawData_3_resample = reSampling_LUT(rawData_3, LUT_procd);
rawData_4_resample = reSampling_LUT(rawData_4, LUT_procd);

% rawData_1   = squeeze(mean(rawDataArray_rescaled(:,6),2));
subplot(2,1,1)
plot(abs(rawData_1(offSet(1):offSet(2)))); hold on
% plot(real(rawData_2));
% plot(real(rawData_3));
plot(real(rawData_4(offSet(1):offSet(2)))); hold off
subplot(2,1,2)
plot(abs(rawData_1_resample(offSet(1):offSet(2)))); hold on
% plot(real(rawData_2));
% plot(real(rawData_3));
plot(real(rawData_4_resample(offSet(1):offSet(2)))); hold off

