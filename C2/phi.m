function phi_val = phi(es)
% This function computes the reduction of resistance factor (phi) of the section
% that, according to ACI318-08, is in function of the strain of the most 
% tensile steel layer (fiber) 
%
% Inputs
% es: strain of the most tensile steel layer
%
% Outputs
% phi: reduction of resistance factor
%
% Notes:
%
    if abs(es) < 0.002
        phi_val = 0.65;
    elseif abs(es) > 0.005
        phi_val = 0.9;
    else
        phi_val = 0.65 + (abs(es)-0.002)*(250/3);
    end
end