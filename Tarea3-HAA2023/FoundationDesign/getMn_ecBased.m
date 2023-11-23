function [M, curvature, c] = getMn_ecBased(Section, N, ec, es_vect)
% Tarea 2 - Hormigón Armado Avanzado
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras R. - Gabriel Ramos V.
%%
% This function determines the nominal flexural strength based on the
% concrete's strain and a range of values for the last reinforcement layer strain
%
% INPUTS
% Section:              Struct of the reinforced concrete's properties
% N:            (kgf)   Axial load
% ec:                   concrete strain in top of section
% es_vect:              vector of es to evaluate
%
% OUTPUTS
% M:            (kgf-cm)
% curvature:    (1/cm)
%
% Notes
%

% Extract variables from Section
fc = Section.fc; % kgf/cm^2
fy = Section.fy; % kgf/cm^2
Es_y = Section.Es;
b = Section.b;
h = Section.h;
% r = frameElement.r;
% nBars = frameElement.nBars;
% diams = frameElement.diams;
% Pu = Section.Pu; % kgf/cm2
% ecu = Section.ecu;
% nLayers = frameElement.nLayers;
% layers = frameElement.layers;
d = Section.d;
as = Section.as;
% P0 = Section.P0;
PC = Section.PC;
beta1_val = Section.beta1_val;

% previous
d_ = max(d);

% get es that produces a P_try almost equal to the axial load N
min_diff = inf;

for i = 1:length(es_vect)   % en vez de ciclo for, puedo ocupar c = ec*d_./(ec+es) y así probar con todos, luego
    % try for a single strain
    es_val = es_vect(i);
    
    % Neutral Axis depth (c)
    c = ec*d_/(ec + es_val);
    
    % Strain of each reinforce layer
    es = ec/c*(c - d); % -
    
    % Sectional Forces
    sigma_Reinf = sigmaReinf(es,fy,Es_y); % kgf/cm^2                          % Steel reinforce bars layers stress
    f_steel = as.*sigma_Reinf; % kgf                                        % Steel reinforce layers forces
    [Cc_vect,  ~] = computeCc(b, h, fc, c, beta1_val); % kgf, cm            % Vector of contributions of each zone of concrete to the nominal flexural strength 
    P_try = sum(Cc_vect) + sum(f_steel);

    % Compare P_try = Cc + sum(f_steel) with N
    if abs(P_try - N) < min_diff
        es_saved = es_val;
        min_diff = P_try;
    end
end
    
% calculate again but with the es that produces the most accurate results
es_val = es_saved;
c = ec*d_/(ec + es_val);
es = ec/c*(c - d); % -
sigma_Reinf = sigmaReinf(es,fy,Es_y); % kgf/cm^2                              % Steel reinforce bars layers stress
f_steel = as.*sigma_Reinf; % kgf                                            % Steel reinforce bars layers forces
[Cc_vect, Cc_centroid] = computeCc(b, h, fc, c, beta1_val); % kgf, cm       % Vector of contributions of each zone of concrete to the nominal flexural strength 
Cc = sum(Cc_vect); % kgf       
    
% results to return
M = Cc*(PC - Cc_centroid) + sum((PC - d).*f_steel); % kgf-cm               % Nominal flexural strength
curvature = ec/c;

end

