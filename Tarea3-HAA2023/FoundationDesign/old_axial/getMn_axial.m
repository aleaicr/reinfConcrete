function [Mn, phi_min] = getMn_axial(Section)
% Tarea 2 - Hormigón Armado Avanzado
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras R. - Gabriel Ramos V.
%%
% This function determines the nominal flexural strength of a section (Mn)
% of a reinforced concrete section
%
% INPUTS
% Section: Struct of the reinforced concrete's properties
% 
% OUTPUTS
% Mn:
% phiMn:
% phi_val:
% ...
%
% Notes
% * Me falta limitar b1c hasta h, porq o si no en es> es_limite el diagrama
% se pasa para el otro lado, pero se puede limitar justo es_limite para que
% de el diagrama sin pasarse e irse a infinito

%% Extract variables from Section
fc = Section.fc; % kgf/cm^2
fy = Section.fy; % kgf/cm^2
Es = Section.Es;
b = Section.b;
h = Section.h;
% r = frameElement.r;
% nBars = frameElement.nBars;
% diams = frameElement.diams;
ecu = Section.ecu;
% nLayers = frameElement.nLayers;
% layers = frameElement.layers;
d = Section.d;
as = Section.as;
% P0 = Section.P0;
PC = Section.PC;
beta1_val = Section.beta1_val;
N = Section.N; % kgf

%% Neutral Axis depth (c)
syms c_
assume(c_, 'integer');
es = (c_ - d)/c_*ecu;
[Cc_vect, ~] = computeCc(b, h, fc, c_, beta1_val); % kgf
eqn = sum(Cc_vect) + sum(as.*sigmaReinf(es,fy,Es)) - N; % kgf, cm
c = double(solve(eqn,c_)); % cm
if ~(c > 0)
    assume(c_ < 0)
    es = (c_ - d)/c_*ecu; % 
    [Cc_vect, ~] = computeCc(b, h, fc, c_, beta1_val); % kgf
    eqn = sum(Cc_vect) + sum(as.*sigmaReinf(es,fy,Es)) - N; % kgf, cm
    c = double(solve(eqn,c_)); % cm
end

%% Reinforcement layers deformation
es = (c - d)/c*ecu; % -

%% Sectional Forces
sigma_Reinf = sigmaReinf(es,fy,Es); % kgf/cm^2                              % Steel reinforce bars layers stress
f_steel = as.*sigma_Reinf; % kgf                                            % Steel reinforce bars layers forces
[Cc_vect, Cc_centroid] = computeCc(b, h, fc, c, beta1_val); % kgf, cm
Cc = sum(Cc_vect); % kgf                                                    % Concrete force

%% Flexural Strength of the Section
Mn = Cc*(PC - Cc_centroid) + sum((PC - d).*f_steel); % kgf-cm               % Nominal flexural strength
phi_min = phi(min(es));                                                     % Strength reduction factor

end

