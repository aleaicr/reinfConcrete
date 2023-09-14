%% Flexure Design Tool
% Tarea 1 -Hormigón Armado Avanzado 2023
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras - Gabriel Ramos

%% Init
clear variables
close all
clc

%% Inputs V25/60 (Viga Preliminar)
% Materials
fc = 300; % kgf/cm2                                                         % Concrete's Strength
fy = 4200; % kgf/cm^2                                                       % Steel Yielding Strength
Es = 2.1*10^6; %kgf/cm^2                                                    % Steel Modulus of Elasticity

% Section geometry
b = 25; % cm                                                                % Width
h = 60; % cm                                                                % Height
r = 5; % cm                                                                 % Recubrimiento

% Reinforcement
nBars = [4; 4];                                                          % Number of bars
diams = [2.2; 2.2]; % cm
% Podemos definid "d" (depths) para la profundidad de cada capa, por si
% vamos a tomar dos tipos de barras para una misma capa (abría que definir
% uno mismo "d" y borrar el que se calcula abajo)

% Axial loading
Pu = 20*1000; % kgf

% ecu
ecu = 0.003;

% lambda
lambda = 1;                                                              % lambda = 1.25 --> Mpr en vez de Mn, lambda = 1 --> Mn

%% Previous Calculations
nLayers = length(diams);                                                    % Number of Layers of longitudinal reinforcement
layers = (1:1:nLayers).';                                                   % Layer IDs
fy = lambda*fy;

% Depth of layers
d = r + (h-2*r)/(nLayers-1)*(layers-1); % cm 

% Area of steel of layers
as = nBars*pi.*(diams/2).^2;    % cm^2

% Axial strength of the section
P0 = 0.85*fc*b*h+sum(as)*(fy - fc); % kgf

% Plastic Centroid
PC = (0.85*b*h^2*fc/2 + sum(as.*d*(fy - fc)))/(P0); % cm

% beta1
beta1_val = beta1(fc);

%% Save Data into Struct
Section = struct();
Section.fc = fc;
Section.fy = fy;
Section.Es = Es;
Section.b = b;
Section.h = h;
Section.r = r;
Section.nBars = nBars;
Section.diams = diams;
Section.Pu = Pu;
Section.ecu = ecu;
Section.nLayers = nLayers;
Section.layers = layers;
Section.d = d;
Section.as = as;
Section.P0 = P0;
Section.PC = PC;
Section.beta1_val = beta1_val;
%% Get Mn
[Mn, phiMn, phi_val, f_steel, Cc, es, c, sigma_Reinf] = getMn(Section); % kgf - cm

%% Display Results
for i = 1:length(nBars)
    fprintf('Refuerzo %.0f: %.0fphi%.0f a %.0f del top\n',i,nBars(i),diams(i)*10,d(i))
end
fprintf('Pu = %.2f [tonf]\n', Pu)
disp('Momento Positivo')
fprintf('\n c = %.2f [cm]\n \n', c)
tabla = table();
tabla.d = d;
tabla.as = as;
tabla.es = es;
tabla.sigma_steel = sigma_Reinf;
tabla.f_steel = f_steel;
disp(tabla)
fprintf('d [cm] | as [cm2]| es [-] | sigma_steel [kgf/cm2] | f_steel [tonf]\n')
fprintf('Cc = %.2f [tonf]\n',Cc/1000)
fprintf('Cc + f_steels- Pu = %.2f [tonf]\n\n', Cc/1000 + sum(f_steel)/1000 - Pu);
tabla = table();
tabla.Mn_tonf_m = Mn/1000/100;
tabla.Phi = phi_val;
tabla.phiMn_tonf_m = phiMn/1000/100;
disp(tabla)

%% Get Mn_neg
[Mn, phiMn, phi_val, f_steel, Cc, es, c, sigma_Reinf] = getMn_neg(Section); % kgf - cm

%% Display Results_neg
disp('Momento Negativo')
fprintf('\n c = %.2f [cm]\n \n', c)
tabla = table();
tabla.d = h-flip(d);
tabla.as = flip(as);
tabla.es = es;
tabla.sigma_steel = sigma_Reinf;
tabla.f_steel = f_steel;
disp(tabla)
fprintf('d [cm] | as [cm2]| es [-] | sigma_steel [kgf/cm2] | f_steel [tonf]\n')
fprintf('Cc = %.2f [tonf]\n',Cc/1000)
fprintf('Cc + f_steels- Pu = %.2f [tonf]\n\n', Cc/1000 + sum(f_steel)/1000 - Pu);
tabla = table();
tabla.Mn_tonf_m = Mn/1000/100;
tabla.Phi = phi_val;
tabla.phiMn_tonf_m = phiMn/1000/100;
disp(tabla)



