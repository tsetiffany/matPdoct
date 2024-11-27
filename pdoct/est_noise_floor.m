%% DOPU processing %%
function [OCT_PN, OCT_SN] = est_noise_floor(cplxData_A,cplxData_B,ref_Frame)

    Nroi = 100;
%     figure(1),
%     subplot(2,1,1),imagesc(imadjust(mat2gray(20*log10(abs(cplxData_A(:,:,ref_Frame)))))),colormap('gray'),title(sprintf('%s Frame %d',patient_id(1:5),ref_Frame));
%     subplot(2,1,2),plot(imadjust(mat2gray(mean(20*log10(abs(cplxData_A(:,:,ref_Frame))),2)))),hold on
%     line([Nroi Nroi],[0 1]),hold off
%     pause(2),
%     saveas(gcf,fullfile(savepath,[patient_id,'_test_refFrame.tif']));
%     close all

    % Noise estimation from reference frame %
    refOCT_P = cplxData_A(:,:,ref_Frame);
    refOCT_S = cplxData_B(:,:,ref_Frame);
    refPhaseOff = repmat(angle(sum(refOCT_S.*conj(refOCT_P),1)), [size(refOCT_S,1) 1]); % new paper; get the value from is/os only? check later
    refOCT_Cplx = refOCT_P + (refOCT_S.*exp(-1j.*refPhaseOff));
    
    % Noise estimate %
    for i = 1:size(refOCT_P,2)
        nOCT_P_real = var(real(refOCT_P(Nroi-10:Nroi+10,i))); % orig 1:20
        nOCT_P_imag = var(imag(refOCT_P(Nroi-10:Nroi+10,i)));
        nOCT_P_cplx(i) = nOCT_P_real + nOCT_P_imag; % Should this be actualy complex noise?????
        nOCT_S_real = var(real(refOCT_S(Nroi-10:Nroi+10,i)));
        nOCT_S_imag = var(imag(refOCT_S(Nroi-10:Nroi+10,i)));
        nOCT_S_cplx(i) = nOCT_S_real + nOCT_S_imag;
    end
    OCT_PN = median(nOCT_P_cplx)+std(nOCT_P_cplx);
    OCT_SN = median(nOCT_S_cplx)+std(nOCT_S_cplx);
%     noise_PS = [OCT_PN, OCT_SN];
%     med_PS=[median(nOCT_P_cplx) median(nOCT_S_cplx)];
%     std_PS=[std(nOCT_P_cplx) std(nOCT_S_cplx)];

end