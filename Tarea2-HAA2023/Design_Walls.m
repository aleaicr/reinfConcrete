%% Diseño Col55/55 Piso 1
% Tarea 2 - Hormigón Armado Avanzado
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras R. - Gabriel Ramos V.

%% Init
clear variables
close all
clc

%% Inputs C55/55 (Columna Piso 1 Exterior)
% Materials
fc = 300; % kgf/cm2                                                         % Concrete's strength
fy = 4200; % kgf/cm^2                                                       % Steel's strength
Es = 2.1*10^6; %kgf/cm^2                                                    % Steel reinforce modulus of elasticity

% Section geometry
b = [300; 50]; % cm                                                         % Concrete section zones widths
h = [50; 300]; % cm                                                         % Concrete section zones heights
% r = 5; % cm

% Reinforcement
diams = [25; 22];   % cm                                                    % Diameter of bars diametro_barras_tipo_1, diam2, diam3
nBars = [6 2; 6 2; 2 0; 2 0; 2 0];                                          % Number of bars per layer
d = [5; 45; 345; 90; 200]; % cm                                             % Depth of each layer

% ecu
ecu = 0.003;

% lambda
lambda = 1;    

% Ultimate Moments and Axial Loads
Mu_ = [5.37; 2.26; 0.97; -5.39; -2.25; -0.94; 0.03; 0.01; 0.07; -0.04; -0.01; -0.05; 0.03; 0.02; 0.08; -0.05; -0.01; -0.04; 5.37; 2.26; 0.97; -5.39; -2.25; -0.94; -0.01; 0.00; 0.02; -0.01; 0.00; 0.02; -0.01; 0.00; 0.02; 4.66; 1.96; 1.11; -4.85; -1.90; -0.80; -0.06; 0.04; 0.21; -0.12; 0.02; 0.10; -0.12; 0.06; 0.31; -0.18; 0.04; 0.21; 4.60; 1.98; 1.21; -4.91; -1.88; -0.69; -0.17; 0.06; 0.29; -0.14; 0.05; 0.24; -0.14; 0.05; 0.23];
Pu_ = -[-265.90; -264.40; -262.90; -274.84; -273.34; -271.84; -237.75; -236.25; -234.76; -302.99; -301.49; -299.99; -399.84; -397.85; -395.85; -465.08; -463.08; -461.08; -427.99; -425.99; -424.00; -436.93; -434.93; -432.94; -476.98; -474.99; -472.99; -420.57; -418.24; -415.92; -400.77; -398.78; -396.78; -238.79; -237.29; -235.79; -263.18; -261.68; -260.18; -215.77; -214.28; -212.78; -286.19; -284.69; -283.20; -365.37; -363.38; -361.38; -435.79; -433.79; -431.80; -388.39; -386.39; -384.39; -412.78; -410.78; -408.78; -441.33; -439.33; -437.34; -390.42; -388.09; -385.76; -371.41; -369.41; -367.41];
Mu2_ = [-1.83; 0.94; 3.75; -2.55; 0.58; 3.68; 1.12; 2.35; 4.09; -5.50; -0.82; 3.34; -0.41; 2.88; 6.69; -7.03; -0.29; 5.94; -3.36; 1.48; 6.35; -4.08; 1.12; 6.28; -4.20; 1.46; 7.13; -3.41; 1.19; 5.78; -3.32; 1.16; 5.64; 2.36; -0.48; -3.11; 1.64; -0.87; -3.60; 5.33; 0.90; -2.98; -1.33; -2.26; -3.73; 6.74; 0.43; -5.34; 0.08; -2.73; -6.09; 3.77; -0.95; -5.46; 3.05; -1.35; -5.96; 3.85; -1.30; -6.45; 3.12; -1.05; -5.22; 3.04; -1.03; -5.09];

% Rango deformaciones 
es_min = -0.00087;                                                          % es mínimo a analizar
es_max = 0.5;                                                               % es máximo a analizar
n_es = 5000;                                                                % Número de puntos dentrod el diagrama de interacciones (notar que no se distribuyen uniformemente)

%% Previous Calculations
% Reinf Layers
nLayers = size(diams);                                                    % Number of Layers of longitudinal reinforcement
layers = (1:1:nLayers).';                                                   % Layer IDs
fy = lambda*fy;

% Concrete zones
nHeight = length(h);
nWidth = length(b);

% Area of concrete for each layer
ag_zones = b.*h; % Area of concrete per zone
ag = sum(ag_zones); % Total area of concrete

% Area of steel for each layer
as_types = pi*(nBars.*(diams.'/10).^2); % cm^2
as = sum(as_types,2);

% Axial Strength of the Section
P0 = 0.85*fc*ag + sum(as)*(fy - fc); % kgf
Ptracc = -sum(as)*fy; % kgf

% Plastic Centroid
PC = (0.85*fc*sum(b.*h.^2)/2 + sum(as.*d*(fy - fc)))/(P0); % cm

% beta1
beta1_val = beta1(fc);

% Fibra más traccionada
d_ = max(d); % solo porq col es simétrica, funciona para ambos lados.

%% Save Data into Struct
Section = struct();
Section.fc = fc;
Section.fy = fy;
Section.Es = Es;
Section.b = b;
Section.h = h;
% Section.r = r;
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
Section.Mu2_ = Mu2_;

%% get Interaction Diagram Data
% Obtener Datos y Graficarlo
[Mn, Pn, phiMn, phiPn] = getInteractionDiagram(Section);
% Obtener datos Mn
Section2 = Section;
Section2.Pu = min(Pu_)*1000;
[Mn_col, phiMn_col, phi_val_col, ~, ~, ~, ~, ~] = getMn(Section2); % kgf y cm
Mn_col = Mn_col/1000/100; % tonf-m
phiMn_col = phiMn_col/1000/100; % tonf-m

%% Display
fprintf('-----------Diseño C55/55-----------\n')
fprintf('-----------Diseño Flexo-Compresión-----------\n')
% Check Acero requerido
Mu = max(max(abs(Mu_)),max(abs(Mu2_)));
As_req = 0.85*fc/fy*b*d_*(1-sqrt(1-2.353*Mu*1000*100/(0.9*b*d_^2*fc)));
fprintf('\nArmadura requerida As_req = %.2f [cm2]\n', As_req)

fprintf('\nArmadura dispuesta-----------\n')
for i = 1:length(nBars)
    fprintf('Refuerzo %.0f: %.0fphi%.0f a %.0f del top, area = %.2f [cm^2]\n',i,nBars(i),diams(i)*10,d(i), as(i))
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
As_min = max(0.8*sqrt(fc)/fy*b*d_,14/fy*b*d_);
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

%check resistencia
fprintf('\nCheck resistencia-----------\n')
if phiMn_col > max(abs(Mu_))
    fprintf('phiMn = %0.2f [tonf] > Mu = %.02f[tonf] OK\n',phiMn_col, max(abs(Mu_)))
else
    fprintf('phiMn = %0.2f [tonf] < Mu = %.02f[tonf] NO OK\n',phiMn_col, max(abs(Mu_)))
end
%% Diseño al corte
fprintf('\n\n-----------DISEÑO CORTE-----------\n\n')

% Inputs
h_col = 500; % cm % Altura de la columna
phi_shear = 0.6; % porque corte viene del análisis
Vu_design = 3.07; % tonf
Vsismico = 0.797; %tonf
fprintf('Vu = %.2f [tonf]\n', Vu_design)
fprintf('Vsismico = %.2f [tonf]\n', Vsismico)

% Vc 
Vc_0 = 0; % si quiero Vc = 0, poner Vc_0=1; si quiero Vc =0.53..., usar Vc_0 = 0
if Vsismico>Vu_design/2 && max(abs(Pu_))<= (b*h)*fc/20
    Vc = 0;
    fprintf('Vc = 0 [tonf]\n')
else
    Vc = 0.53*sqrt(fc)*b*d_/1000;
    fprintf('Vc = %.2f [tonf]\n',Vc)
end

if Vc_0 == 1
    fprintf('Vamos a ocupar Vc = 0[tonf]\n')
    Vc = 0;
end

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
s = 12;
phi_tr = 8; % mm
hx = (max(h,b)-2*r)/2; %% Cambiar la formula o el valor, esta formula es solo porq funciona en esta configuración de estribos + trabas
Av_s_tabla = 4.19; %cm2/m
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
