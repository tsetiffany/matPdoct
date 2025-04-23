function attcoef = attenuation_coefficient(im)
% Calculate attenuation coefficient using Vermeer2014BOE model.
% axial_pix_len = 1;
axial_pix_len = 5.26 / 1000;

weight_fac = cumsum(im,1);

compfac = weight_fac - im;
compfac(compfac==0) = nan;

intensity_ratio = im ./ compfac;
intensity_ratio(intensity_ratio<=0)=0;
intensity_ratio(isnan(intensity_ratio))=0;

attcoef = log(1+intensity_ratio) / (2*axial_pix_len);

% figure(1),
% subplot(1,2,1),imagesc(imadjust(mat2gray(im(51:500,:)))), colormap('gray')
% subplot(1,2,2),imagesc(imadjust(mat2gray(attcoef(51:500,:)))), colormap('gray')

end