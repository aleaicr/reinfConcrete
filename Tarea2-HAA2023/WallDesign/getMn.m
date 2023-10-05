function [Mn, phiMn, phi_val, f_steel, Cc, es, c, sigma_Reinf] = getMn(Section)
% This function determines the nominal flexural strength of a section (Mn)
% of a reinforced concrete section
%
% INPUTS
% Section: Struct of the reinforced concrete's properties
% 
% OUTPUTS
%

% Extract variables from Section
fc = Section.fc; % kgf/cm^2
fy = Section.fy; % kgf/cm^2
Es = Section.Es;
b = Section.b;
% h = Section.h;
% r = frameElement.r;
% nBars = frameElement.nBars;
% diams = frameElement.diams;
Pu = Section.Pu;
ecu = Section.ecu;
% nLayers = frameElement.nLayers;
% layers = frameElement.layers;
d = Section.d;
as = Section.as;
% P0 = Section.P0;
PC = Section.PC;
beta1_val = Section.beta1_val;

% Neutral Axis depth (c)
syms c_
assume(c_, 'positive');
es = (c_ - d)/c_*ecu;
[Cc_vect, ~] = computeCc(b, h, fc, c_, beta1_val);
eqn = sum(Cc_vect) + sum(as.*sigmaReinf(es,fy,Es)) - Pu;
c = double(solve(eqn,c_)); % cm
if ~(c > 0)
    assume(c_, 'negative')
    es = (c_ - d)/c_*ecu;
    [Cc_vect, ~] = computeCc(b, h, fc, c_, beta1_val);
    eqn = sum(Cc_vect) + sum(as.*sigmaReinf(es,fy,Es)) - Pu;
    c = double(solve(eqn,c_)); % cm
end

% Reinforcement layers deformation
es = (c - d)/c*ecu; % -

% Sectional Forces
sigma_Reinf = sigmaReinf(es,fy,Es); % kgf/cm^2                              % Steel reinforce bars layers stress
f_steel = as.*sigma_Reinf; % kgf                                            % Steel reinforce bars layers forces
[Cc_vect, Cc_centroid] = computeCc(b, h, fc, c, beta1_val); % kgf, cm
Cc = sum(Cc_vect); % kgf                                                    % Concrete force

% Flexural Strength of the Section
Mn = Cc*(PC - Cc_centroid) + sum((PC - d).*f_steel); % kgf-cm             % Nominal flexural strength
phi_val = phi(min(es));                                                     % Strength reduction factor
phiMn = phi_val*Mn; % kgf-cm                                                % Design flexural strength

end

