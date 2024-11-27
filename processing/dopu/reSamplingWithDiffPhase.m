function rescaledData = reSamplingWithDiffPhase(fftData, rescaleParam, diffPhase, offSetIdx)

fftData_R = padarray(real(fftData), [(2^15 - size(fftData,1)) 0], 'post');
fftData_I = padarray(imag(fftData), [(2^15 - size(fftData,1)) 0], 'post');
fftData_pad = fftData_R + (1j.*fftData_I);



rawData_real        = real(ifft(fftData_pad));
rawData_imag        = imag(ifft(fftData_pad));
rawData_rescaled    = zeros([length(rescaleParam), size(fftData,2)]);
rawData_corrected   = zeros([length(rescaleParam), size(fftData,2)]);
idxPixel            = linspace(1,size(fftData,1),size(fftData,1));


for idxData = 1:size(fftData,2)
    resPoints = rescaleParam * (size(fftData_pad,1)/size(fftData,1));
    rawData_rescaled = (interp1([1:size(rawData_real, 1)]',rawData_real(:,idxData),resPoints,'spline'))...
        +1j.*(interp1([1:size(rawData_imag, 1)]',rawData_imag(:,idxData),resPoints,'spline'));
    fftData_rescaled = fft(rawData_rescaled);
    [~, maxIdx]        = sort(abs(fftData_rescaled(offSetIdx:offSetIdx+100)), 'descend');
    fftData_corrected  = fftData_rescaled.*exp(-1j*(idxPixel'*(diffPhase(idxData)/(maxIdx(1)+offSetIdx-1))));
    rawData_corrected(:,idxData) = (ifft(fftData_corrected));

end

rescaledData = rawData_corrected;

end


