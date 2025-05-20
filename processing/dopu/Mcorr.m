        usfac = 200;
        for FrameNum = 1 : numBscans
            ProcdData_AB(: , : , FrameNum) = abs(ProcdData_A(: , : , FrameNum)) + abs(ProcdData_B(: , : , FrameNum));
        end
        ref_OCT_A = imgaussfilt(abs(ProcdData_A(:, :, 20)), 2);
        % ref_OCT_B = imgaussfilt(abs(volume_mcorr_ChP(:, :, 20)), 2);
        [volume_mcorr_PS, volume_mcorr_ChP, volume_mcorr_ChS] = globalReg_PS(ProcdData_AB , ProcdData_A , ProcdData_B , usfac , ref_OCT_A);
