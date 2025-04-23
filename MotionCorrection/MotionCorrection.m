function vol_mcorr = MotionCorrection(vol)

    %% Global axial mcorr
    usfac = 20; % upsampling factor
    
    % For iterative mcorr
    vol_mcorr = vol;
    
%     figure;
%     subplot(2,2,1); imagesc(mean(abs(vol_mcorr),3)); title('Unaligned averaged Bscan'); xlabel("A-scans"); ylabel("Depth"); set(gca, FontSize=15)
    
    % Iterate axial correction until max error is below 0.25 pixels
    ii = 0;
    yshift_global = 1;
    while max(abs(yshift_global)) > 0.25 % Iterates until maximum axial shift is < 0.25 px
        ii = ii + 1;
        [vol_mcorr, yshift_global, xshift_global] = mcorrLocal_axial(vol_mcorr, usfac);
        if ii > 10
            break
        end
    end
%     disp("Axial registration completed in " + num2str(ii) + " iterations")
    
    % Plot the results
%     subplot(2,2,2); imagesc(mean(abs(vol_mcorr),3)); title('Aligned averaged Bscan'); xlabel("A-scans"); ylabel("Depth"); set(gca, FontSize=15)
%     subplot(2,2,3); plot(yshift_global); title('Axial (y) pixel shifts'); xlabel("Bscans"); ylabel("y shift magnitude"); set(gca, FontSize=15)
%     subplot(2,2,4); plot(xshift_global); title('Lateral (x) pixel shifts'); xlabel("Bscans"); ylabel("x shift magnitude"); set(gca, FontSize=15)
    
    % Lateral correction with evolving reference frame
    [vol_mcorr, xshift] = mcorrLocal_lateral(vol_mcorr, usfac);
    
    % Plot the lateral correction curve
%     figure, plot(xshift)

end