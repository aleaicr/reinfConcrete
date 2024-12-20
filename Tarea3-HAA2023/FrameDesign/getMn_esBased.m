function [Mn, Pn, phi_min] = getMn_esBased(Section)
% Tarea 2 - Hormigón Armado Avanzado
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras R. - Gabriel Ramos V.
%%
% This function determines the Mn of a reinforced concrete section based on
% the desired "most deformed steel layer" parameter 'es'
%
% Inputs
%
% Outputs
%
% Notes
%
%%
% Extract variables from Section
fc = Section.fc; % kgf/cm^2
fy = Section.fy; % kgf/cm^2
Es = Section.Es;
b = Section.b;
h = Section.h;
% r = frameElement.r;
% nBars = frameElement.nBars;
% diams = frameElement.diams;
% Pu = Section.Pu;
ecu = Section.ecu;
% nLayers = frameElement.nLayers;
% layers = frameElement.layers;
d = Section.d;
as = Section.as;
% P0 = Section.P0;
PC = Section.PC;
beta1_val = Section.beta1_val;
es_val = Section.es_val;
d_ = max(d);

% Neutral Axis depth (c)
c = ecu*d_/(ecu + es_val);

% Reinforcement layers strain
es = ecu/c*(c - d); % -

% Sectional Forces
sigma_Reinf = sigmaReinf(es,fy,Es); % kgf/cm^2                              % Steel reinforce bars layers stress
f_steel = as.*sigma_Reinf; % kgf                                            % Steel reinforce bars layers forces
[Cc_vect, Cc_centroid] = computeCc(b, h, fc, c, beta1_val); % kgf, cm       % Vector of contributions of each zone of concrete to the nominal flexural strength 
Cc = sum(Cc_vect); % kgf                                                    % Concrete force

% Flexural Strength of the Section
Mn = Cc*(PC - Cc_centroid) + sum((PC - d).*f_steel); % kgf-cm               % Nominal flexural strength
Pn = Cc + sum(f_steel);
phi_min = phi(min(es));                                                     % Strength reduction factor
% phi_curvature = ecu/c;                                                      % Curvature

end

