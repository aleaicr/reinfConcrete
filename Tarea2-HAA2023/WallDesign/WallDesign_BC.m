%% Wall Design (T-shaped wall)
% Tarea 2 - Hormigón Armado Avanzado
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras R. - Gabriel Ramos V.
%
% Notes
% * The functions getInteractionDiagram, getMn, getMn_esBased; are evaluated
% for this T-shaped wall only (I never tried more sections), but I think
% that all is in function of the vectors/matrices sizes in the input
%% Init
clear variables
close all
clc

%% Inputs T-shaped Wall
% Materials
fc = 300; % kgf/cm2                                                         % Concrete's strength
fy = 4200; % kgf/cm^2                                                       % Steel's strength
Es = 2.1*10^6; %kgf/cm^2                                                    % Steel reinforce modulus of elasticity

% Section geometry
b = [790; 50]; % cm                                                         % Concrete section zones widths
h = [30; 640]; % cm                                                         % Concrete section zones heights

% Reinforcement
diams = [20; 12];   % cm                                                     % Diameter of bars diametro_barras_tipo_1, diam2, diam3
nBars = [2 17; 2 17; 0 2; 0 2; 0 2; 0 2; 0 2; 0 2; 0 2; 0 2; 0 2; 0 2; 0 2; 0 2; 0 2; 0 2; 0 2; 0 2; 0 2];                                                         % Number of bars per layer
d = [35; 70; 105; 140; 175; 210; 245; 280; 315; 350; 385; 420; 455; 490; 525; 560; 595; 630; 665];

% ecu
ecu = 0.003;

% lambda
lambda = 1;    % lambda*fy (puede tomar 1 o 1.25 si quiero calcular Mpr)    % No confundir con lambda de la ACI318 el cual corresponde a factor de reducción por hormigón "ligero" (lightweight)

% Pier internal forces Data (from etabs analysis results tables) (tonf, m)
% col1: P, col2: V2, col3: V3, col4: T, col5: M2, col6: M3
% each row is a load combination
% for wall design, we're only interested in the first floor internal forces
%
internalForcesData = [
-1473.7327 348.0111 145.0541 646.2986 3502.5691 4092.9903;
-1953.787 -340.0257 -157.3547 -647.0617 -3452.2185 -4097.0131;
-1575.7248 7.2082 340.7293 4.4981 5961.4339 2.5527;
-1851.7949 0.7772 -353.0299 -5.2612 -5911.0832 -6.5755;
-2523.1071 9.9282 336.6747 4.1342 5977.0811 1.2399;
-2799.1772 3.4971 -357.0845 -5.6252 -5895.4361 -7.8883;
-2421.115 350.7311 140.9995 645.9346 3518.2163 4091.6775;
-2901.1693 -337.3057 -161.4094 -647.4257 -3436.5713 -4098.3258;
-2893.9367 7.5493 -11.4085 -0.8832 45.4478 -3.714;
-2665.8487 6.2109 -9.5671 -0.5935 39.1616 -3.1288;
-1473.7327 350.7311 340.7293 646.2986 5977.0811 4092.9903;
-2901.1693 -340.0257 -357.0845 -647.4257 -5911.0832 -4098.3258;
-2495.8525 6.0283 -9.2053 -0.6131 38.0654 -3.0172
];

Mu_ = internalForcesData(:,6);                                              % As we're interested in M3 now
Pu_ = -internalForcesData(:,1);

% Deformation range
% for interaction diagram
es_min = -0.0002;                                                           % es mínimo a analizar
es_max = 0.5;                                                               % es máximo a analizar
n_es = 5000;                                                                % Número de puntos dentrod el diagrama de interacciones (notar que no se distribuyen uniformemente)

% For moment-curvature diagram
ec_min = 0.000001;
ec_max = 0.003;
n_ec = 50;

% for interaction diagram (axial based)
n_N = 200; % interpolation points between pure compression and pure traction

% concrete partitions
part = 200;

%% Previous Calculations
% Reinf Layers
nLayers = size(diams);                                                      % Number of Layers of longitudinal reinforcement
layers = (1:1:nLayers).';                                                   % Layer IDs
fy = lambda*fy;

% Concrete zones
nHeight = length(h);            % Notar que nHeight = nWidth = nZones
nWidth = length(b);

% Area of concrete for each layer
ag_zones = b.*h; % cm^2                                                     % Area of concrete per zone
ag = sum(ag_zones); % cm^2                                                  % Total area of concrete

% Area of steel for each layer
as_types = pi*(nBars.*(diams.'/10).^2); % cm^2
as = sum(as_types,2); % cm^2                                                % Vector of as in each layer of reinforcement

% Axial Strength of the Section
P0 = 0.85*fc*ag + sum(as)*(fy - fc); % kgf
Ptracc = -sum(as)*fy; % kgf

% Plastic Centroid
PC = (0.85*fc*sum(b.*h.^2)/2 + sum(as.*d*(fy - fc)))/(P0); % cm

% beta1
beta1_val = beta1(fc);

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
Section.ecc = [ec_min; ec_max; n_ec];
Section.Mu_ = Mu_;
Section.Pu_ = Pu_;
Section.n_N = n_N;
Section.part = part;

%% get Interaction Diagram Data
% graficar diagrama de interacción y momento-curvatura
[Mn, Pn, phiMn, phiPn] = getInteractionDiagram(Section);
% [Mn, Pn, phiMn, phiPn] = getInteractionDiagram_axial(Section);
% [M, curvature, M_neg, curvature_neg] = getMomentCurvature(Section);


%% Diseño a flexión
fprintf('-----------------Diseño a flexión-----------------')
fprintf('Cuantía------------------------------------------\n')



% %% Display
% fprintf('-----------Wall Design-----------\n')
% fprintf('Note that this Wall is a T-Shaped Wall\n\n')
% fprintf('-----------Flexure and Axial Design-----------\n')
% fprintf('\nArmadura dispuesta-----------\n')
% for i = 1:length(nBars)
%     fprintf('Refuerzo %.0f: %.0fphi%.0f + %.0fphi%.0f a %.0f del top, area = %.2f [cm^2]\n',i, nBars(i,1), diams(1),nBars(i,2), diams(2), d(i), as(i))
% end
% % Check Geomería
% 
% % Check Rango Cuantia
% fprintf('\ncheck cuantia_min > 0.0025-----------\n')
% cuantia = sum(as)/(b*h);
% if cuantia > 0.0025
%     fprintf('Cuantia = %.4f OK\n', cuantia);
% else
%     fprintf('Cuantia = %.2f NO OK\n', cuantia);
% end
% 
% % check hx
% fprintf('\nCheck espaciamiento barras longitudinales-----------\n')
% esp_long = max(diff(d));  % como es simétrica aplica para ambos sentidos
% if esp_long < 15
%     fprintf('esp_long = %0.2f [cm] < 15[cm] OK\n',esp_long)
% else
%     fprintf('esp_long = %0.2f [cm] > 15[cm] NO OK\n',esp_long)
% end
% 
% %check resistencia
% fprintf('\nCheck resistencia-----------\n')
% if phiMn_col > max(abs(Mu_))
%     fprintf('phiMn = %0.2f [tonf] > Mu = %.02f[tonf] OK\n',phiMn_col, max(abs(Mu_)))
% else
%     fprintf('phiMn = %0.2f [tonf] < Mu = %.02f[tonf] NO OK\n',phiMn_col, max(abs(Mu_)))
% end
% %% Diseño al corte
% fprintf('\n\n-----------DISEÑO CORTE-----------\n\n')
% 
% % Inputs
% h_col = 500; % cm % Altura de la columna
% phi_shear = 0.6; % porque corte viene del análisis
% Vu_design = 3.07; % tonf
% Vsismico = 0.797; %tonf
% fprintf('Vu = %.2f [tonf]\n', Vu_design)
% fprintf('Vsismico = %.2f [tonf]\n', Vsismico)
% 
% % Vc 
% Vc_0 = 0; % si quiero Vc = 0, poner Vc_0=1; si quiero Vc =0.53..., usar Vc_0 = 0
% if Vsismico>Vu_design/2 && max(abs(Pu_))<= (b*h)*fc/20
%     Vc = 0;
%     fprintf('Vc = 0 [tonf]\n')
% else
%     Vc = 0.53*sqrt(fc)*b*d_/1000;
%     fprintf('Vc = %.2f [tonf]\n',Vc)
% end
% 
% if Vc_0 == 1
%     fprintf('Vamos a ocupar Vc = 0[tonf]\n')
%     Vc = 0;
% end
% 
% % Steel
% Vs_req = Vu_design/phi_shear - Vc; % tonf
% Av_s_req = Vs_req*1000/(fy*d_); % cm^2/cm
% Av_s_min = max(0.25*sqrt(fc)*b/fy,3.5*b/fy);
% fprintf('Vs_req = %.2f [tonf]\n',Vs_req)
% fprintf('Av_s_req = %.3f [cm2/cm]\n', Av_s_req)
% fprintf('Av_s_min = %.3f [cm2/cm]\n', Av_s_min)
% 
% %%%% CONFIGURACIÓN ARMADURA Armadura
% estribos = 1;
% trabas = 1;
% nbarras = 2*estribos + 1*trabas;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Av_s_buscar
% Av_s_buscar = max(Av_s_req*100/nbarras,Av_s_min*100/nbarras);
% fprintf('\nBuscar Av_s = %.2f [cm2/m]\n',Av_s_buscar)
% fprintf('Buscar Av_s = %.4f [cm2/cm]\n\n',Av_s_buscar/100)
% %%%%%% ELEGIR ACÁ LA ARMADURA, LUEGO DE BUSCAR
% 
% phi_estr = 8; %mm
% s = 12;
% phi_tr = 8; % mm
% hx = (max(h,b)-2*r)/2; %% Cambiar la formula o el valor, esta formula es solo porq funciona en esta configuración de estribos + trabas
% Av_s_tabla = 4.19; %cm2/m
% % NO modificar acá abajo
% s_estr = s; % cm
% s_tr = s; %cm
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Av_s = nbarras*Av_s_tabla/100; % cm2/cm
% Av = estribos*2*pi*0.25*(phi_estr/10)^2 + trabas*pi*0.25*(phi_tr/10)^2;
% s = max(s_tr,s_estr);
% Vs = fy*d_*Av_s/1000; % tonf
% 
% % display dispuesto
% fprintf('Av_s_dispuesto = %.3f [cm2/m]\n',Av_s*100)
% fprintf('Av_s_dispuesto = %.3f [cm2/cm]\n',Av_s)
% fprintf('s_dispuesto = %.2f [cm]\n',s)
% fprintf('Vs_dispuesto = %.2f [tonf]\n',Vs)
% 
% % Check Resistencia
% fprintf('\nCheck Resistencia-----------\n')
% Vn = Vs + Vc; % tonf
% phiVn = phi_shear*Vn;
% if phiVn > Vu_design
%     fprintf('phiVn = %.2f [Tonf] > Vu = %.2f [tonf] OK \n',phiVn,Vu_design)
% else
%     fprintf('phiVn = %.2f [Tonf] < Vu = %.2f [tonf] NO OK \n',phiVn,Vu_design)
% end
% 
% % Check cuantía
% fprintf('\nCheck cuantía Av_s-----------\n')
% if Av_s > Av_s_min && Av_s > Av_s_req
%     fprintf('Av_s > Av_s_min y req OK\n')
% else
%     fprintf('Cuantía NO OK \n')
% end
% 
% % Check Vsmax
% fprintf('\nCheck Vs_max-----------\n')
% Vs_max = 2.2*sqrt(fc)*b*d_/1000;
% if Vs < Vs_max
%     fprintf('Vs < Vs_max = %.2f[tonf] OK\n',Vs_max)
% else
%     fprintf('Vs > Vs_max = %.2f[tonf] NO OK\n',Vs_max)
% end
% 
% % S_max
% fprintf('\nCheck Espaciamiento s_max-----------\n')
% % 11.4.5
% fprintf('11.4.5 - Spacing limits for shear reinforcement\n')
% Vs_limit_s = 1.1*sqrt(fc)*b*d_/1000; % tonf
% fprintf('Vs_limit = 1.1sqrt(fc)bd = %.2f [tonf]\n',Vs_limit_s)
% if Vs < Vs_limit_s % tonf
%     s_max_1145 = d_/2; % cm
% else
%     s_max_1145 = d_/4; % cm
% end
% s0 = min(max(10 + (35-hx)/3, 10), 15); % cm
% s_max_21643 = min([max(h,b)/4, min(h,b)/4, 6*min(diams), s0]); % cm
% s_max_21352 = min([max(h,b)/2, min(h,b)/2, 8*min(diams), 30]); % cm
% s_max_710 = min([16*min(diams), 48*min(phi_estr,phi_tr)/10, max(h,b), min(h,b)]); % cm
% s_max_21645 = min(6*min(diams),15); % cm
% s1 = min([s_max_21643; s_max_21352; s_max_710; s_max_1145]); % cm
% s2 = min([s_max_21645: s_max_710, s_max_1145]); %cm
% l0 = max([h; b; h_col/6; 45]);
% fprintf('l0 = %.1f [cm]\n',l0)
% fprintf('s1 = %.3f [cm]\n',s1)
% fprintf('s2 = %.3f [cm]\n',s2)
% 
% if s < s2
%     fprintf('s = %.2f [cm] < s2 = %.2f [cm]  OK\n',s,s2)
% else
%     fprintf('s = %.2f [cm] > s2 = %.2f [cm]  NO OK\n',s,s2)
% end
% if s < s1
%     fprintf('s = %.2f [cm] < s1 = %.2f [cm]  OK \n',s,s1)
% else
%     fprintf('s = %.2f [cm] > s1 = %.2f [cm]  NO OK\n',s,s1)
% end
% 
% % A disponer
% s1_disp = min(s, redondearEsp(s1));
% s2_disp = min(s, redondearEsp(s2));
% [l0_dispuesto,nl0] = l0_dispuestoFunc(l0, s1_disp);
% fprintf('Se dispondrán E%.0fphi%.0f@%.0f + TRphi%.0f@%.0f fuera de l0\n',estribos, phi_estr,s2_disp,phi_tr,s2_disp)
% fprintf('Dentro de l0=%.0f [cm] se dispondrán %.0f refuerzos E%.0fphi%.0f@%.0f + TRphi%.0f@%.0f\n', l0_dispuesto, nl0, estribos, phi_estr, s1_disp, phi_tr, s1_disp)
