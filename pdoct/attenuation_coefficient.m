function attcoef = attenuation_coefficient(im)
% Calculate attenuation coefficient using Vermeer2014BOE model.
% axial_pix_len = 1;
axial_pix_len = 5.26 / 1000;

weight_fac = cumsum(im,1);

compfac = weight_fac - im;
compfac(compfac==0) = nan;

intensity_ratio = im ./ compfac;
intensity_ratio(intensity_ratio<=0) = 0;
intensity_ratio(isnan(intensity_ratio)) = 0;

attcoef = log(1 + intensity_ratio) / (2 * axial_pix_len);

end
