%% Diseño Col55/55 Piso 1
% Tarea 1 - Hormigón Armado Avanzado
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras R. - Gabriel Ramos V.

%% Init
clear variables
close all
clc

%% Inputs C55/55 (Columna Piso 1 Exterior)
% Materials
fc = 300; % kgf/cm2
fy = 4200; % kgf/cm^2
Es = 2.1*10^6; %kgf/cm^2

% Section geometry
b = 65; % cm
h = 65; % cm
r = 5; % cm

% Reinforcement
nBars = [6; 2; 2; 2; 2; 6];                                                             % Number of bars per layer
diams = [3.2; 3.2; 3.2; 3.2; 3.2; 3.2]; % cm                                                    % Diameter of bars in the layer
% Podemos definir "d" (depths) para la profundidad de cada capa, por si
% vamos a tomar dos tipos de barras para una misma capa (abría que definir
% uno mismo "d" y borrar el que se calcula abajo)

% ecu
ecu = 0.003;

% lambda
lambda = 1;    

% Ultimate Moments and Axial Loads
Mu_ = [8.81; 4.48; 1.40; -9.48; -4.24; -0.27; 0.01; 0.24; 1.16; -0.68; 0.00; -0.02; -0.21; 0.32; 1.53; -0.90; 0.07; 0.36; 8.59; 4.55; 1.78; -9.70; -4.16; 0.11; -0.62; 0.22; 1.06; -0.52; 0.18; 0.88; -0.50; 0.18; 0.86; 7.81; 4.47; 7.39; -12.78; -2.72; 1.08; -2.15; 0.99; 4.80; -2.81; 0.76; 3.67; -3.86; 1.59; 7.71; -4.52; 1.36; 6.58; 6.11; 5.07; 10.30; -14.48; -2.11; 3.99; -4.72; 1.66; 8.05; -3.86; 1.36; 6.59; -3.76; 1.33; 6.41];
Pu_ = -[-429.56; -427.47; -425.38; -456.49; -454.40; -452.31; -427.59; -425.50; -423.41; -458.46; -456.37; -454.28; -706.81; -704.02; -701.23; -737.68; -734.89; -732.10; -708.78; -705.99; -703.20; -735.71; -732.92; -730.13; -803.38; -800.59; -797.80; -689.15; -685.90; -682.65; -663.55; -660.76; -657.97; -231.92; -229.83; -227.74; -425.54; -423.45; -421.36; -296.53; -294.44; -292.35; -360.93; -358.83; -356.74; -500.86; -498.07; -495.29; -565.25; -562.47; -559.68; -436.25; -433.46; -430.67; -629.87; -627.08; -624.29; -591.40; -588.61; -585.82; -511.36; -508.10; -504.85; -490.44; -487.65; -484.87];
Mu2_ = [1.37; 0.15; -0.69; -0.16; -0.53; -1.29; 7.04; 2.89; -0.25; -5.83; -3.27; -1.73; 7.45; 2.76; -0.92; -5.41; -3.41; -2.41; 1.78; 0.02; -1.37; 0.26; -0.66; -1.96; 1.15; -0.36; -1.87; 0.94; -0.30; -1.54; 0.91; -0.29; -1.49; 0.40; 0.48; 0.56; -0.97; -0.23; 0.49; 6.13; 3.21; 1.24; -6.69; -2.97; -0.19; 5.93; 3.30; 1.61; -6.89; -2.88; 0.18; 0.21; 0.56; 0.93; -1.16; -0.15; 0.85; -0.54; 0.23; 1.01; -0.44; 0.19; 0.82; -0.42; 0.18; 0.79];

% Rango deformaciones 
es_min = -0.00087;
es_max = 0.5;
n_es = 5000;

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
P0 = 0.85*fc*b*h + sum(as)*(fy - fc); % kgf
Ptracc = -sum(as)*fy; % kgf

% Plastic Centroid
PC = (0.85*b*h^2*fc/2 + sum(as.*d*(fy - fc)))/(P0); % cm

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
Section.Mu2_ = Mu2_;

%% get Interaction Diagram Data
% Obtener Datos y Graficarlo
[Mn, Pn, phiMn, phiPn] = getInteractionDiagram(Section);
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
Vu_design = 5.63; % tonf
Vsismico = 3.05; %tonf
fprintf('Vu = %.2f [tonf]\n', Vu_design)
fprintf('Vsismico = %.2f [tonf]\n', Vsismico)

% Vc 
if Vsismico>Vu_design/2 && max(abs(Pu_))<= (b*h)*fc/20
    Vc = 0;
    fprintf('Vc = 0 [tonf]\n')
else
    Vc = 0.53*sqrt(fc)*b*d_/1000;
    fprintf('Vc = %.2f [tonf]\n',Vc)
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
s = 14;
phi_tr = 8; % mm
hx = (max(h,b)-2*r)/5*3; %% Cambiar la formula o el valor, esta formula es solo porq funciona en esta configuración de estribos + trabas
Av_s_tabla = 3.59; %cm2/m
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
