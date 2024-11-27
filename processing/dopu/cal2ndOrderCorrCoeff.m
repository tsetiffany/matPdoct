function [numer, denom, denom_norm] = cal2ndOrderCorrCoeff(hermitX_R, hermitX_I, Int_A, Int_B, noise, SNR_Correct)

hermitX_cplx = hermitX_R + (i.*hermitX_I);
numer        = abs(hermitX_cplx)./noise;

if SNR_Correct == 1
    array_A     = Int_A - noise;
    array_B     = Int_B - noise;    
    array_C     = (array_A - array_B).^2;
    array_D     = (array_A + array_B).^2;
    array_E     = 1 - (array_C./array_D);
    array_E(array_E>1) = 1;
    array_E(array_E<0) = 0;
 
    denom       = (sqrt(array_E).*((array_A + array_B)/2))./noise;
    denom_norm  = sqrt(Int_A.*Int_B)./noise;

else
    
    denom       = sqrt(Int_A.*Int_B)./noise;
    denom_norm  = denom;

end