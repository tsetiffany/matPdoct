function [vol_mcorr, yshift, xshift] = mcorrLocal_axial(vol, usfac)

    numPoints   = size(vol, 1);
    numLines    = size(vol, 2);
%     numFrames   = size(vol, 3);   % For global motion correction
    numFrames   = 3;       % For local motion correction
    
    kLinear         = linspace(-1,1,numPoints);
    kaxis           = repmat(kLinear',1,numLines);
    
    vol_mcorr = vol;
    
    reference = abs(vol(:, :, 1));
    for I = 2:numFrames
        %%% For every Bscan, reference frame will be the previous Bscan %%%

        [output, Greg, cplx_shift, diffphase] = dftregistration_AA(fft2(abs(reference)),...
            fft2(abs(vol(:,:,I))), usfac);
        yshift(I) = output(3);
        xshift(I) = output(4);

        % This line does the lateral shift
%         temp = ifft2(fft2(vol(:,:,I)).*cplx_shift);

        %% Only do axial correction
        vol_mcorr(:,:,I) = fft(ifft(vol_mcorr(:,:,I)).*exp(2j.*(output(3)*(kaxis))));

    end
end