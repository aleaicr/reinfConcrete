function beta1_val = beta1(fc)
% This functions compute the value of Beta1 of the ACI318-08
%
% INPUTS
% fc: [kgf/cm2] concrete's strength
%
% OUTPUTS
% beta1_val: reduction factor of the height of the compressed concrete zone
%
% Notes
    fcMPa = fc/10;          % fc in MPa = fc in kgf/cm2 / 10 according to ACI318-08
    if fcMPa >= 17 && fcMPa <= 28
        beta1_val = 0.85;
    elseif fcMPa > 28 && fcMPa < 55
        beta1_val = 0.85 - (0.05/7)*(fcMPa-28);
    else
        beta1_val = 0.65;
    end
end
