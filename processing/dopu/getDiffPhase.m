function phaseSlope = getDiffPhase(fftData_Cal, rescaleParam, diffPhase)

fftData_R = padarray(real(fftData_Cal), [(2^15 - size(fftData_Cal,1)) 0], 'post');
fftData_I = padarray(imag(fftData_Cal), [(2^15 - size(fftData_Cal,1)) 0], 'post');
fftData_pad = fftData_R + (1j.*fftData_I);

rawData_real        = real(ifft(fftData_pad));
rawData_imag        = imag(ifft(fftData_pad));
rawData_rescaled    = zeros([length(rescaleParam), size(fftData_Cal,2)]);    
idxPixel            = linspace(1,size(fftData_Cal,1),size(fftData_Cal,1));


for idxData = 1:size(fftData_Cal,2)
    resPoints = rescaleParam * (size(fftData_pad,1)/size(fftData_Cal,1));

    rawData_rescaled(:,idxData) = (interp1([1:size(rawData_real, 1)]',rawData_real(:,idxData),resPoints,'spline'))...
        +1j.*(interp1([1:size(rawData_imag, 1)]',rawData_imag(:,idxData),resPoints,'spline'));

    fftData_rescaled(:,idxData) = fft(rawData_rescaled(:,idxData));
    
    [~, maxIdx]        = sort(abs(fftData_rescaled(:, idxData)), 'descend');
    phaseSlope(:,idxData) = (idxPixel'*(diffPhase(idxData)/maxIdx(1)));

end

end


