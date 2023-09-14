%% Interaction Diagram for a Section
% Tarea 1 - Hormigón Armado Avanzado
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras R. - Gabriel Ramos V.

%% Init
clear variables
close all
clc

%% Inputs V25/60 (Viga Preliminar)
% Materials
fc = 300; % kgf/cm2
fy = 4200; % kgf/cm^2
Es = 2.1*10^6; %kgf/cm^2

% Section geometry
b = 30; % cm
h = 55; % cm
r = 5; % cm

% Reinforcement
nBars = [4; 4];                                                             % Number of bars per layer
diams = [2.2; 2.2]; % cm                                                    % Diameter of bars in the layer
% Podemos definid "d" (depths) para la profundidad de cada capa, por si
% vamos a tomar dos tipos de barras para una misma capa (abría que definir
% uno mismo "d" y borrar el que se calcula abajo)

% ecu
ecu = 0.003;

% lambda
lambda = 1;    

% Ultimate Moments and Axial Loads
Mu_ = [12; 5; 15; 5; 15; 20]; 
Pu_ = [20; 30; 15; 20; 50; 50];

es_min = -0.00088;
es_max = 0.5;
n_es = 1000;

%% Previous Calculations
% Layers
nLayers = length(diams);                                                    % Number of Layers of longitudinal reinforcement
layers = (1:1:nLayers).';                                                   % Layer IDs
fy = lambda*fy;

% Depth of each layer
d = r + (h-2*r)/(nLayers-1)*(layers-1); % cm 

% Area of steel for each layer
as = nBars*pi.*(diams/2).^2;    % cm^2

% Axial Strength of the Section
P0 = 0.85*fc*b*h+sum(as)*(fy - fc); % kgf
Ptracc = -sum(as)*fy; % kgf

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
% Section.Pu = Pu;
Section.ecu = ecu;
Section.nLayers = nLayers;
Section.layers = layers;
Section.d = d;
Section.as = as;
Section.P0 = P0;
Section.PC = PC;
Section.beta1_val = beta1_val;
Section.ess = [es_min; es_max; n_es];
Section.Mu_ = Mu_;
Section.Pu_ = Pu_;

%% get Interaction Diagram Data
[Mn, Pn, phiMn, phiPn] = getInteractionDiagram(Section);

% Notar que getInteractionDiagram dibuja además el diagrama de interacción,
% no es necesario volver a dibujarlo acá

