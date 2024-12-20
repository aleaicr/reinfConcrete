function [Mn, Pn, phi_val, phi_curvature] = getMn_esBased(Section)
% This function determines the Mn of a reinforced concrete section based on
% the desired "most deformed steel layer" parameter 'es'

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
es_eval = Section.es_eval;

% Neutral Axis depth (c)
c = max(d)*ecu/(ecu+es_eval);

% Reinforcement layers deformation
es = (c - d)/c*ecu; % -

% Sectional Forces
sigma_Reinf = sigmaReinf(es,fy,Es); % kgf/cm^2                              % Steel reinforce bars layers stress
f_steel = as.*sigma_Reinf; % kgf                                            % Steel reinforce bars layers forces
Cc = 0.85*fc*beta1_val*c*b; % kgf                                           % Concrete force

% Flexural Strength of the Section
Mn = Cc*(PC - beta1_val*c/2) + sum((PC - d).*f_steel); % kgf-cm             % Nominal flexural strength
Pn = Cc + sum(f_steel);
phi_val = phi(min(es));                                                     % Strength reduction factor
phi_curvature = ecu/c;

end

