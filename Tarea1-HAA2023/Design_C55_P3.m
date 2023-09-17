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
b = 55; % cm
h = 55; % cm
r = 5; % cm

% Reinforcement
nBars = [5; 2; 2; 2; 5];                                                             % Number of bars per layer
diams = [2.2; 2.2; 2.2; 2.2; 2.2]; % cm                                                    % Diameter of bars in the layer
% Podemos definir "d" (depths) para la profundidad de cada capa, por si
% vamos a tomar dos tipos de barras para una misma capa (abría que definir
% uno mismo "d" y borrar el que se calcula abajo)

% ecu
ecu = 0.003;

% lambda
lambda = 1;    

% Ultimate Moments and Axial Loads
Mu_ = [7.9371; 1.9232; 4.4889; -8.2557; -1.9514; -4.2267; 0.059; 0.0095; 0.3022; -0.3775; -0.0377; -0.04; -0.0388; 0.0011; 0.3831; -0.4753; -0.0461; 0.0409; 7.8394; 1.9148; 4.5698; -8.3535; -1.9598; -4.1458; -0.2858; -0.025; 0.2359; -0.2478; -0.0219; 0.2039; -0.2411; -0.0213; 0.1984; 6.50; 1.53; 4.14; -7.27; -1.63; -3.56; -0.20; -0.03; 0.43; -0.57; -0.07; 0.15; -0.46; -0.06; 0.63; -0.83; -0.10; 0.34; 6.24; 1.50; 4.33; -7.53; -1.66; -3.37; -0.72; -0.09; 0.54; -0.60; -0.07; 0.45; -0.58; -0.07; 0.44];
Pu_ = -[-222.3558; -221.3178; -220.2799; -230.5342; -229.4962; -228.4583; -196.6974; -195.6595; -194.6215; -256.1925; -255.1546; -254.1166; -331.7874; -330.4034; -329.0195; -391.2825; -389.8985; -388.5146; -357.4457; -356.0618; -354.6779; -365.6241; -364.2402; -362.8562; -398.6406; -397.2567; -395.8727; -352.2478; -350.6332; -349.0186; -336.0211; -334.6371; -333.2532; -196.93; -195.89; -194.86; -219.65; -218.61; -217.58; -176.14; -175.10; -174.06; -240.45; -239.41; -238.37; -299.56; -298.18; -296.79; -363.87; -362.49; -361.10; -320.36; -318.97; -317.59; -343.08; -341.69; -340.31; -365.30; -363.92; -362.53; -324.01; -322.39; -320.78; -308.52; -307.14; -305.75];
Mu2_ = [-5.9763; -0.8117; 4.8039; -7.0795; -1.1388; 4.351; -2.1611; 0.3244; 6.4537; -10.8947; -2.2748; 2.7012; -6.7554; -0.3593; 9.6806; -15.489; -2.9586; 5.928; -10.5706; -1.4954; 8.0307; -11.6738; -1.8225; 7.5779; -12.5756; -1.8744; 8.8268; -10.1545; -1.517; 7.1205; -9.9207; -1.4805; 6.9598; 6.40; 1.03; -3.18; 4.66; 0.70; -4.42; 10.03; 2.18; -1.83; 1.02; -0.45; -5.78; 13.96; 2.79; -4.54; 4.95; 0.16; -8.49; 10.33; 1.64; -5.90; 8.59; 1.31; -7.13; 10.71; 1.66; -7.38; 8.60; 1.34; -5.92; 8.41; 1.31; -5.79];

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
h_col = 365; % cm % Altura de la columna
phi_shear = 0.6; % porque corte viene del análisis
Vu_design = 8.24; % tonf
Vsismico = 2.03; %tonf
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
s = 12;
phi_tr = 8; % mm
hx = (max(h,b)-2*r)/2; %% Cambiar la formula o el valor, esta formula es solo porq funciona en esta configuración de estribos + trabas
Av_s_tabla = 4.19; %cm2/m
% 8,14 = 3.59; 8,12= 4.59; 8, 10=5.03
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
