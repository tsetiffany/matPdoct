%% Reference Frame process %%
% ref_frame = round(numBscans/2);
function [dispCoeffs, ref_FFTData_HamWin,ref_FFT_DisComp]  = est_dispersion_coeff(fn,LUT,parameters,ref_frame,dispROI,dispMaxOrder)
%{
Input
-----
fn:
    file name
LUT:
    resampline parameters
parameters:
    acqusiton parameters

Output
-----
dispCoeffs:
    dispersion coefficients
%}
    %%% Preset parameter %%%
    coeffRange      = 50;
    bitDepth        = 12;
    byteSize        = bitDepth/8;

    %%% Acquisition Parameters %%%
%     parameters    = getParameters(fn);
    numPoints     = parameters(1)*2;
    numAscans     = parameters(2)/2;
    numBscans     = parameters(3);
    numCscans     = parameters(4);
    numMscans     = parameters(5);

    %%% Load reference frame %%%
%     ref_frame = round(numBscans/2)+1;
    if bitDepth == 16
        fid = fopen(fn);
        fseek(fid,byteSize*numPoints*numAscans*(ref_frame-1),-1);
        ref_RawData_interlace = fread(fid,[numPoints,numAscans], 'uint16');
        fclose(fid);
    elseif bitDepth == 12
        ref_RawData_interlace = unpack_u12u16(fn,numPoints,numAscans,ref_frame);
    else
        error('%d-bit data not supported.', bitDepth)
    end
    ref_RawData = hilbert(cat(2, ref_RawData_interlace(1:2:end,:),ref_RawData_interlace(2:2:end,:)));
    %figure(1),plot(real(ref_RawData(:,1))), hold on, plot(real(ref_RawData(:,end/2+1))),hold off,
    
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
    %figure(1),imagesc(imadjust(mat2gray(20.*log10(abs(ref_FFTData_FPNSub(51:750,:)))))); colormap(gray)
    
    %%% Dispersion Estimation %%%
%     dispROI        = [100, 500];

    dispCoeffs = setDispCoeff(ref_RawData_HamWin,dispROI,dispMaxOrder,coeffRange);
    ref_RawData_DisComp = compDisPhase(ref_RawData_HamWin,dispMaxOrder,dispCoeffs);
    ref_FFT_DisComp   = fft(ref_RawData_DisComp);
%     ref_FFT_DisComp     = 20.*log10(abs(ref_FFT_Final));
    %figure(2),imagesc(imadjust(mat2gray(ref_OCT_Log(51:750,:)))); colormap(gray)

end
