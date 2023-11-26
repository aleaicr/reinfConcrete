function beta1_val = beta1(fc)
    fcMPa = fc/10;
    if fcMPa >= 17 && fcMPa <= 28
        beta1_val = 0.85;
    elseif fcMPa > 28 && fcMPa < 55
        beta1_val = 0.85 - (0.05/7)*(fcMPa-28);
    else
        beta1_val = 0.65;
    end
end
