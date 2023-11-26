%% Wall Design (T-shaped wall)
% Tarea 2 - Hormigón Armado Avanzado
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras R. - Gabriel Ramos V.
%
% Notes
% * The functions getInteractionDiagram, getMn, getMn_esBased; are evaluated
% for this T-shaped wall only (I never tried more sections), but I think
% that all is in function of the vectors/matrices sizes in the input
%

%% Init
clear variables
close all
clc

%% Inputs T-shaped Wall
% Materials
fc = 300; % kgf/cm2                                                         % Concrete's strength
fy = 4200; % kgf/cm^2                                                       % Steel's strength
Es = 2.1*10^6; %kgf/cm^2                                                    % Steel reinforce modulus of elasticity

% lambda
lambda = 1;    % lambda*fy (puede tomar 1 o 1.25 si quiero calcular Mpr)    % No confundir con lambda de la ACI318 el cual corresponde a factor de reducción por hormigón "ligero" (lightweight)

% Section geometry
b = 70; % cm                                                                % Concrete section zones widths
h = 70; % cm                                                                % Concrete section zones heights
r = 5;

% Reinforcement (from top to bottom) ITERATION 2
diams = [8; 8; 8; 8; 8; 8; 8];   % mm                                       % diameters of each type of bar
nBars = [7; 2; 2; 2; 2; 2; 7];                                              % number of bars of each type in each layer

% ecu
ecu = 0.003;

% Strain range for interaction diagram
es_min = -0.0005;                                                           % es mínimo a analizar
es_max = 1;                                                                 % es máximo a analizar
n_es = 50000;                                                               % Número de puntos dentrod el diagrama de interacciones (notar que no se distribuyen uniformemente)

% Columns IDs
% columnsIDs = [17; 32; 18; 31; 19; 30]; % 65x65
columnsIDs = [15; 16; 13; 14]; % 70x70

% Dir. of the SAP200 results file
fileDir = '../ModeloconVF.xlsx';

%% Previous calculations
% Reinf Layers
nLayers = length(diams);                                                    % Number of Layers of longitudinal reinforcement
layers = (1:1:nLayers).';                                                   % Layer IDs
fy = lambda*fy;

% Depth of each layer
d = r + (h-2*r)/(nLayers-1)*(layers-1); % cm 

% Concrete zones
nHeight = length(h);            % Notar que nHeight = nWidth = nZones
nWidth = length(b);

% Area of concrete for each layer
ag_zones = b.*h; % cm^2                                                     % Area of concrete per zone
ag = sum(ag_zones); % cm^2                                                  % Total area of concrete

% Area of steel for each layer
as_types = 0.25*pi*(nBars.*(diams.'/10).^2); % cm^2
as = sum(as_types,2); % cm^2                                                % Vector of as in each layer of reinforcement

% Axial Strength of the Section
P0 = 0.85*fc*(ag-sum(as)) + sum(as*fy); % kgf
Ptracc = -sum(as)*fy; % kgf

% Plastic Centroid
aux = h;
aux(end) = [];
aux = [0; cumsum(aux)];
h_centroids = h/2 + aux; % cm
PC = (0.85*fc*sum(b.*h.*h_centroids) + sum(as.*d*(fy-0.85*fc)))/P0; % cm

% beta1
beta1_val = beta1(fc);

%% Loads
% get all internal forces for LRFD cases for each column
% loadComb = 'LRFD';
% for i = 1: length(columnsIDs)
%     [internalLoads, allTable] = getFrameLoads(fileDir, columnsIDs, loadComb);
% end
% 
% % All columns internal loads in one matrix
% intLoads = zeros(1,6);
% for i = 1:length(internalLoads)
%     intLoads = [intLoads; internalLoads(i).frameTable{:,7:12}];             % Can't find another way to store all matrices
% end
% intLoads(1, :) = [];                                                        % Delete first row
intLoads = readmatrix('C70x70Loads.csv');                                % Ya lo hice, dejé comentado porq se demora mucho en leer todas
Mu_ = [intLoads(:, 6); intLoads(:, 5)]; % tonf-m
Pu_ = -[intLoads(:,1); intLoads(:,1)]; % tonf-m
Vu_ = [intLoads(:,3); intLoads(:,2)]; % tonf-m

%% Save Data into Struct
Section = struct();
Section.fc = fc;
Section.fy = fy;
Section.Es = Es;
Section.b = b;
Section.h = h;
% Section.r = r;
% Section.nBars = nBars;
% Section.diams = diams;
% Section.Pu = Pu;
Section.ecu = ecu;
% Section.nLayers = nLayers;
% Section.layers = layers;
Section.d = d;
Section.as = as;
Section.P0 = P0;
Section.PC = PC;
Section.beta1_val = beta1_val;
Section.ess = [es_min; es_max; n_es];
% Section.ess_2 = [es_min_2; es_max_2; n_es_2];
% Section.ecc = [ec_min; ec_max; n_ec];
Section.Mu_ = Mu_;
Section.Pu_ = Pu_;
% Section.n_N = n_N;
% Section.part = part;

%% get Diagrams
[Mn, Pn, phiMn, phiPn] = getInteractionDiagram(Section);

%% DISEÑO FRAME
%% Display
fprintf('-----------Diseño C70/70-----------\n')
fprintf('-----------Diseño Flexo-Compresión-----------\n')
% Check Acero requerido
% Mu = max(max(abs(Mu_)),max(abs(Mu2_)));
Mu = max(abs(Mu_));
d_ = max(d);
As_req = 0.85*fc/fy*b*d_*(1-sqrt(1-2.353*Mu*1000*100/(0.9*b*d_^2*fc)));
fprintf('\nArmadura requerida As_req = %.2f [cm2]\n', As_req)

fprintf('\nArmadura dispuesta-----------\n')
for i = 1:length(nBars)
    fprintf('Refuerzo %.0f: %.0fphi%.0f a %.0f del top, area = %.2f [cm^2]\n',i,nBars(i),diams(i),d(i), as(i))
end
% Check Geomería
fprintf('\nCheck geometría-----------\n')
if min(b,h) > 30
    fprintf('b > 30cm OK\n')
else
    fprintf('b < 30cm NO OK\n')
end
if min(b,h)/max(b,h) > 0.
    fprintf('b/h > 0.4 OK\n')
else
    fprintf('b/h < 0.4 NO OK\n')
end

% Check Rango Cuantia
fprintf('\ncheck cuantia entre 0.01 y 0.06-----------\n')
cuantia = sum(as)/(b*h);
if cuantia > 0.01 && cuantia < 0.06
    fprintf('Cuantia = %.4f OK\n', cuantia);
else
    fprintf('Cuantia = %.2f NO OK\n', cuantia);
end

% Cuantia minima
fprintf('\ncheck cuantia mayor a minima-----------\n')
As_min = min([0.005; max(0.8*sqrt(fc)/fy*b*d_,14/fy*b*d_)]);
cuantia_min = As_min/(b*h);
if cuantia > cuantia_min
    fprintf('Cuantia mayor a cuantia minima = %0.4f OK\n',cuantia_min)
else
    fprintf('Cuantia no es mayor a cuantía minima= %0.4f no OK\n', cuantia_min)
end

% check hx
fprintf('\nCheck espaciamiento barras longitudinales-----------\n')
esp_long = max(diff(d));  % como es simétrica aplica para ambos sentidos
if esp_long < 15
    fprintf('esp_long = %0.2f [cm] < 15[cm] OK\n',esp_long)
else
    fprintf('esp_long = %0.2f [cm] > 15[cm] NO OK\n',esp_long)
end

% %check resistencia --> Diagrama de interacción

%% Diseño al corte
fprintf('\n\n-----------DISEÑO CORTE-----------\n\n')

% Inputs
h_col = 130; % cm % Altura de la columna
phi_shear = 0.6; % porque corte viene del análisis
Vu_design = max(abs(Vu_)); % tonf
fprintf('Vu = %.2f [tonf]\n', Vu_design)

% Vc
Vc = 0.53*sqrt(fc)*b*d_/1000;
fprintf('Vc = %.2f [tonf]\n',Vc)

% Steel
Vs_req = Vu_design/phi_shear - Vc; % tonf
Av_s_req = Vs_req*1000/(fy*d_); % cm^2/cm
Av_s_min = max(0.25*sqrt(fc)*b/fy,3.5*b/fy);
fprintf('Vs_req = %.2f [tonf]\n',Vs_req)
fprintf('Av_s_req = %.3f [cm2/cm]\n', Av_s_req)
fprintf('Av_s_min = %.3f [cm2/cm]\n', Av_s_min)

%%%% CONFIGURACIÓN ARMADURA Armadura
estribos = 1;
trabas = 1;
nbarras = 2*estribos + 1*trabas;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Av_s_buscar
Av_s_buscar = max(Av_s_req*100/nbarras,Av_s_min*100/nbarras);
fprintf('\nBuscar Av_s = %.2f [cm2/m]\n',Av_s_buscar)
fprintf('Buscar Av_s = %.4f [cm2/cm]\n\n',Av_s_buscar/100)
%%%%%% ELEGIR ACÁ LA ARMADURA, LUEGO DE BUSCAR

phi_estr = 8; %mm
s = 10;
phi_tr = 8; % mm
hx = (max(h,b)-2*r)/2; %% Cambiar la formula o el valor, esta formula es solo porq funciona en esta configuración de estribos + trabas
Av_s_tabla = 5.03; %cm2/m
% 8@14 = 3.59; 8@12= 4.19; 8@10=5.03, 8@20=2.51
% NO modificar acá abajo
s_estr = s; % cm
s_tr = s; %cm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Av_s = nbarras*Av_s_tabla/100; % cm2/cm
Av = estribos*2*pi*0.25*(phi_estr/10)^2 + trabas*pi*0.25*(phi_tr/10)^2;
s = max(s_tr,s_estr);
Vs = fy*d_*Av_s/1000; % tonf

% display dispuesto
fprintf('Av_s_dispuesto = %.3f [cm2/m]\n',Av_s*100)
fprintf('Av_s_dispuesto = %.3f [cm2/cm]\n',Av_s)
fprintf('s_dispuesto = %.2f [cm]\n',s)
fprintf('Vs_dispuesto = %.2f [tonf]\n',Vs)

% Check Resistencia
fprintf('\nCheck Resistencia-----------\n')
Vn = Vs + Vc; % tonf
phiVn = phi_shear*Vn;
if phiVn > Vu_design
    fprintf('phiVn = %.2f [Tonf] > Vu = %.2f [tonf] OK \n',phiVn,Vu_design)
else
    fprintf('phiVn = %.2f [Tonf] < Vu = %.2f [tonf] NO OK \n',phiVn,Vu_design)
end

% Check cuantía
fprintf('\nCheck cuantía Av_s-----------\n')
if Av_s > Av_s_min && Av_s > Av_s_req
    fprintf('Av_s > Av_s_min y req OK\n')
else
    fprintf('Cuantía NO OK \n')
end

% Check cuantía mínima
rho_min_1 = 0.09*fc/fy;
rho_min_2 = 0.3*fc/fy*(ag/(h-r)^2-1);
rho_min = min([rho_min_1;rho_min_2]);
fprintf('Avs_min_21.6.4.4 = %.4f\n y rho_disp = %.4f',rho_min, Av_s)

% Check Vsmax
fprintf('\nCheck Vs_max-----------\n')
Vs_max = 2.2*sqrt(fc)*b*d_/1000;
if Vs < Vs_max
    fprintf('Vs < Vs_max = %.2f[tonf] OK\n',Vs_max)
else
    fprintf('Vs > Vs_max = %.2f[tonf] NO OK\n',Vs_max)
end

% S_max
fprintf('\nCheck Espaciamiento s_max-----------\n')
% 11.4.5
fprintf('11.4.5 - Spacing limits for shear reinforcement\n')
Vs_limit_s = 1.1*sqrt(fc)*b*d_/1000; % tonf
fprintf('Vs_limit = 1.1sqrt(fc)bd = %.2f [tonf]\n',Vs_limit_s)
if Vs < Vs_limit_s % tonf
    s_max_1145 = d_/2; % cm
else
    s_max_1145 = d_/4; % cm
end
s0 = min(max(10 + (35-hx)/3, 10), 15); % cm
s_max_21643 = min([max(h,b)/4, min(h,b)/4, 6*min(diams), s0]); % cm
s_max_21352 = min([max(h,b)/2, min(h,b)/2, 8*min(diams), 30]); % cm
s_max_710 = min([16*min(diams), 48*min(phi_estr,phi_tr)/10, max(h,b), min(h,b)]); % cm
s_max_21645 = min(6*min(diams),15); % cm
s1 = min([s_max_21643; s_max_21352; s_max_710; s_max_1145]); % cm
s2 = min([s_max_21645: s_max_710, s_max_1145]); %cm
l0 = max([h; b; h_col/6; 45]);
fprintf('l0 = %.1f [cm]\n',l0)
fprintf('s1 = %.3f [cm]\n',s1)
fprintf('s2 = %.3f [cm]\n',s2)

if s < s2
    fprintf('s = %.2f [cm] < s2 = %.2f [cm]  OK\n',s,s2)
else
    fprintf('s = %.2f [cm] > s2 = %.2f [cm]  NO OK\n',s,s2)
end
if s < s1
    fprintf('s = %.2f [cm] < s1 = %.2f [cm]  OK \n',s,s1)
else
    fprintf('s = %.2f [cm] > s1 = %.2f [cm]  NO OK\n',s,s1)
end

% A disponer
s1_disp = min(s, redondearEsp(s1));
s2_disp = min(s, redondearEsp(s2));
[l0_dispuesto,nl0] = l0_dispuestoFunc(l0, s1_disp);
fprintf('Se dispondrán E%.0fphi%.0f@%.0f + TRphi%.0f@%.0f fuera de l0\n',estribos, phi_estr,s2_disp,phi_tr,s2_disp)
fprintf('Dentro de l0=%.0f [cm] se dispondrán %.0f refuerzos E%.0fphi%.0f@%.0f + TRphi%.0f@%.0f\n', l0_dispuesto, nl0, estribos, phi_estr, s1_disp, phi_tr, s1_disp)
