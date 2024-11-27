function diffPhases = estPhaseDiff(refAline, fftData, offSetIdx)

ref_fftData_R   = padarray(real(refAline), [(2^15 - size(refAline,1)) 0], 'post');
ref_fftData_I   = padarray(imag(refAline), [(2^15 - size(refAline,1)) 0], 'post');
ref_fftData_pad = ref_fftData_R + (1j.*ref_fftData_I);

tar_fftData_R   = padarray(real(fftData), [(2^15 - size(fftData,1)) 0], 'post');
tar_fftData_I   = padarray(imag(fftData), [(2^15 - size(fftData,1)) 0], 'post');
tar_fftData_pad = tar_fftData_R + (1j.*tar_fftData_I);

refCalSig_raw = ifft(ref_fftData_pad);
tarCalSigs_raw = ifft(tar_fftData_pad);

refCalSig      = fft(refCalSig_raw);
tarCalSigs     = fft(tarCalSigs_raw);

cplxConjCalSig  = tarCalSigs.*repmat(conj(refCalSig), [1 size(tarCalSigs,2)]);
[~, idx]        = sort(abs(refCalSig(offSetIdx:end/2)), 'descend');
refMaxIndx      = offSetIdx+idx(1)-1;
diffPhases      = angle(cplxConjCalSig(refMaxIndx,:));


end


