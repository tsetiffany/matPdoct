function [vol_mcorr, xshift] = mcorrLocal_lateral(vol, usfac)

    numPoints   = size(vol, 1);
    numLines    = size(vol, 2);
%     numFrames   = size(vol, 3);     % For global motion correction
    numFrames   = 3;       % For local motion correction
    
    kLinear         = linspace(-1,1,numPoints);
    kaxis           = repmat(kLinear',1,numLines);

    xshift = zeros([numFrames 1]);
    xshiftlimit = 20;
    count = 0;
    ref_reset = 0;
    count_reset = 1;
    reset_ind = 0;
    
    vol_mcorr = vol;
    
    % Start with reference as 1st Bscan
    reference = abs(vol_mcorr(:, :, 1));

    I = 2;
    while I <= numFrames
        %%% For every Bscan, reference frame will be the previous Bscan %%%

        [output, Greg, cplx_shift, diffphase] = dftregistration_AA(fft2(abs(reference)),...
            fft2(abs(vol_mcorr(:,:,I))), usfac);
        xshift(I) = output(4);

        % If the X shifts are outside limits, reset the reference image 
%         if abs(xshift(I)) >= xshiftlimit || abs(xshift(I)-xshift(I-1)) > 1
        if abs(xshift(I)-xshift(I-1)) > 1  
            ref_reset = ref_reset + 1;
            reference = vol_mcorr(:, :, I-1);
            count = count + 1;
            reset_ind(count_reset) = I;
            count_reset = count_reset + 1;
            disp(count_reset)
            if count > 1
                % If reseting the reference didn't work, apply the previous
                % shift to this Bscan, then use this Bscan as the new
                % reference

                % Use the shift from previous Bscan
                xshift(I) = xshift(I-1);

                % Use the relevant part of dftregistration with predefined shifts
                [nr,nc]=size(fft2(abs(vol_mcorr(:,:,I))));
                Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
                Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
                [Nc,Nr] = meshgrid(Nc,Nr);

                % Apply the lateral shift
                cplx_shift = exp(1i*2*pi*(-xshift(I)*Nc/nc)); 
                vol_mcorr(:,:,I) = ifft2(fft2(vol_mcorr(:,:,I)).*cplx_shift);

                I = I + 1; % Move on to the next Bscan
                count = 0;
                pause(0.01);
            end
            continue

            I = I; % Do not move on to the next Bscan yet
        end
        count = 0;

        %% Apply the lateral shift
        vol_mcorr(:,:,I) = ifft2(fft2(vol_mcorr(:,:,I)).*cplx_shift);

        I = I + 1; % Move on to the next Bscan

    end
end