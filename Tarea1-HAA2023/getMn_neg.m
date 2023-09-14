function [Mn, phiMn, phi_val, f_steel, Cc, es, c, sigma_Reinf] = getMn_neg(Section)
% This function determines the Mn of a reinforced concrete section

% Extract variables from frameData
fc = Section.fc; % kgf/cm^2
fy = Section.fy; % kgf/cm^2
Es = Section.Es;
b = Section.b;
h = Section.h;
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

% Flip section
d = h-flip(d);
as = flip(as);
PC = h-PC;

% Neutral Axis depth (c)
syms c_
assume(c_, 'positive');
es = (c_ - d)/c_*ecu;
eqn = 0.85*fc*beta1_val*c_*b + sum(as.*sigmaReinf(es,fy,Es)) - Pu;
c = double(solve(eqn,c_)); % cm
if ~(c > 0)
    assume(c_, 'negative')
    es = (c_ - d)/c_*ecu;
    eqn = 0.85*fc*beta1_val*c_*b + sum(as.*sigmaReinf(es,fy,Es)) - Pu;
    c = double(solve(eqn,c_)); % cm
end

% Reinforcement layers deformation
es = (c - d)/c*ecu; % -

% Sectional Forces
sigma_Reinf = sigmaReinf(es,fy,Es); % kgf/cm^2                              % Steel reinforce bars layers stress
f_steel = as.*sigma_Reinf; % kgf                                            % Steel reinforce bars layers forces
Cc = 0.85*fc*beta1_val*c*b; % kgf                                           % Concrete force

% Flexural Strength of the Section
Mn = Cc*(PC - beta1_val*c/2) + sum((PC - d).*f_steel); % kgf-cm             % Nominal flexural strength
phi_val = phi(min(es));                                                     % Strength reduction factor
phiMn = phi_val*Mn; % kgf-cm                                                % Design flexural strength

end