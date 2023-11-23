%% Foundation Flexural Design
% Tarea 3 - Hormigón Armado Avanzado
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

% lambda
lambda = 1;    % lambda*fy (puede tomar 1 o 1.25 si quiero calcular Mpr)    % No confundir con lambda de la ACI318 el cual corresponde a factor de reducción por hormigón "ligero" (lightweight)

% Section geometry
b = 100; % cm                                                               % Concrete section zones widths
h = 70; % cm                                                                % Concrete section zones heights
r = 7; % cm                                                                 % Reinforcement cover
s = 10; % cm                                                                % Spacing

% Reinforcement
diams = 25;  % mm                                                            % Reinforcement diameter
nBars = (b - s)/s + 1;                                                          % Number of bars
d = h - r; % cm                                                             % Effective depth

% ecu
ecu = 0.003;

% % Pier internal forces Data
% % Copy paste etabs analysis results tables from excel into a .txt file
% % col1: P, col2: V2, col3: V3, col4: T, col5: M2, col6: M3
% internalForcesFileName = 'resultsEtabsPiers.xlsm'; % data must be in tonf, m
% Piers = {'PA', 'PB', 'PC', 'PD'};
% stories = [1; 2; 3];
% 
% % Strain range for interaction diagram
% es_min = -0.0005;                                                           % es mínimo a analizar
% es_max = 1;                                                                 % es máximo a analizar
% n_es = 50000;                                                                % Número de puntos dentrod el diagrama de interacciones (notar que no se distribuyen uniformemente)
% 
% % Strain range for moment-curvature diagram
% ec_min = 0.00001;
% ec_max = 0.003;
% n_ec = 20;
% es_min_2 = -2;  % intentar no modificar
% es_max_2 = 2;   % intentar no modificar
% n_es_2 = 50000; % aumentar si no da suficientemente definido el M-C
% 
% % Number of axial loads to make the moment-curvature diagram, will be
% % equally spaced within min(Pu_) and max(Pu_)
% N_partitions = 2;  
% 
% % Demand of displacement (at the top of the building)
% du = 30; % cm

% % concrete partitions for use in mandel model (not available yet)
% part = 200;

%% Previous Calculations
% load internal forces
% [~, ~, ~, internalForcesData] = readAnalysisResults(internalForcesFileName, Piers, stories);
% clc
% Mu_ = internalForcesData(:,6);  % M3   % tonf-m                             % As we're interested in M3 now
% Pu_ = -internalForcesData(:,1); % P    % ton
% Pu_ = zeros(length(Mu_),1); % tonf                                           % Neglecting the axial load in the flexural design
Pu_ = 0;
Mu_ = 0;

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

%% Save section data into a struct
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
% Section.ess = [es_min; es_max; n_es];
% Section.ess_2 = [es_min_2; es_max_2; n_es_2];
% Section.ecc = [ec_min; ec_max; n_ec];
Section.Mu_ = Mu_;
Section.Pu_ = Pu_;
% Section.n_N = n_N;
% Section.part = part;

%% get Mn
Section.N = 0;
[Mn, phiMn, phi_val, f_steel, Cc, es, c, sigma_Reinf, phi_curvature] = getMn(Section); % kgf - cm (inputs y otuputs)

% transform to tonf - m
Mn = Mn/1000/100; % tonf-m/m                                              % Porque b = 100cm = 1m
phiMn = phiMn/1000/100; % tonf-m/m

fprintf('Mn = %.3f [tonf-m/m]\n', Mn)
fprintf('phiMn = %.3f [tonf-m/m]\n', phiMn)

% f_steel = f_steel/1000; % tonf
% Cc = Cc/1000; % tonf
% sigma_Reinf = sigma_Reinf/1000/100^2;

%% Wall Flexural Design
% % plot Interaction Diagram
% [Mn, Pn, phiMn, phiPn] = getInteractionDiagram(Section);
% Mn, phiMn, tonf-
% Pn, phiPn, tonf
% % plot Moment Curvature Diagram
% Solo temporalmente
% ecc_old = Section.ecc;
% Section.ecc = [ec_min; ec_max; 5];
% [M, curvature, c, M_neg, curvature_neg, c_neg] = getMomentCurvature(Section, N_partitions);
% Section.ecc = ecc_old;
% M, M_neg, tonf-m
% curvature, curvature_neg, 1/cm
% c, c_neg, cm
% 
% %% Diseño a flexión
% fprintf('-----------------Diseño a flexión-----------------\n')
% fprintf('Cuantía------------------------------------------\n')
% rho_l = sum(as)/ag; % -
% rho_l_min = 1.8/1000; % -
% if rho_l > rho_l_min
%     fprintf('Cuantía rho_l = %.4f > rho_l_min = %.4f OK\n',rho_l,rho_l_min)
% else
%     fprintf('Cuantía rho_l = %.4f < rho_l_min = %.4f NO OK\n',rho_l,rho_l_min)
% end
% fprintf('S_max------------------------------------------\n')
% s = max(diff(d)); % cm
% s_max = 45; % cm
% if s < s_max
%     fprintf('s = %.1f < s_max = %.1f OK\n',s,s_max)
% else
%     fprintf('s = %.1f > s_max = %.1f NO OK\n',s,s_max)
% end
% fprintf('Pu_max------------------------------------------\n')
% Pu_lim = 0.35*ag*fc/1000; % tonf
% Pu_max = max(Pu_); % tonf
% if Pu_max < Pu_lim
%     fprintf('Pu_max = %.0f < Pu_lim = 0.35Agf = %.1f OK\n',Pu_max,Pu_lim)
% else
%     fprintf('Pu_max = %.0f > Pu_lim = 0.35Agfc = %.1f NO OK\n',Pu_max,Pu_lim)
% end
% 
% %% Demanda de curvatura
% fprintf('-----------------Demanda de curvatura-----------------\n')
% fprintf('Curvatura mínima al 0.008------------------------------------------\n')
% h_w = 4515; % cm    % distance from the critical section to the top of the building
% l_w = sum(h); % cm
% phi_u = 2*du/(h_w*l_w)*100; % 1/m
% fprintf('phi_u = 2du/(hwlw) = %.6f [1/m] \n',phi_u)
% c_pos = max(c(:,end)); %  cm
% c_neg_ = max(abs(c_neg(:,end))); % cm
% ec_demanda_pos = phi_u/100*c_pos;
% ec_demanda_neg = phi_u/100*c_neg_;
% 
% if ec_demanda_pos < 0.008
%     fprintf('ec_demandado_pos = phi_u*c = %.4f < 0.008 OK\n', ec_demanda_pos)
% else
%     fprintf('ec_demandado_pos = phi_u*c = %.4f > 0.008 NO OK \n', ec_demanda_pos)
%     fprintf('por lo tanto se requiere elemento de borde "abajo"\n')
% end
% 
% if ec_demanda_neg < 0.008
%     fprintf('ec_demandado_neg = phi_u*c = %.4f < 0.008 OK\n', ec_demanda_neg)
% else
%     fprintf('ec_demandado_neg = phi_u*c = %.4f > 0.008 NO OK \n', ec_demanda_neg)
%     fprintf('por lo tanto se requiere elemento de borde "arriba"\n')
% end
% 
% % En el gráfico que se va a generar, se ocupa ec_último = 0.008, si la 
% % curva NO supera las líneas (horizontalmente), entonces tiene que cambiar 
% % la disposición del refuerzo, si la % supera, entonces se podría necesitar
% % elemento de borde en función si el siguiente gráfico (siguiente
% % comentario). Si no lo supera, entonces, hay que agregarle, si lo supera,
% % entonces no hay que hacer nada
% ec_min_2 = 0.00001;
% ec_max_2 = 0.008;
% n_ec_2 = n_ec;
% Section2 = Section;
% Section2.ecc = [ec_min_2; ec_max_2; n_ec_2];
% [M_2, curvature_2, c_2, M_neg_2, curvature_neg_2, c_neg_2] = getMomentCurvature(Section2, 1);
% hold on
% xline(phi_u,'--r', '$$\phi_u$$', 'Interpreter', 'latex', 'fontsize', 20,'linewidth', 3)
% xline(-phi_u,'--r', '$$-\phi_u$$', 'Interpreter', 'latex', 'fontsize', 20, 'linewidth', 3)
% hold off
% legend(sprintf('P_u = %.1f [tonf]', max(Pu_)));
% 
% 
% % En el gráfico que se va a generar, se ocupa ec_último = 0.003, si la
% % curva NO supera los \phi_u (lineas verticales), entonces se necesita un
% % elemento de borde en ese sentido para cubrir lo que falta para alcanzar
% % ec = 0.008 de deformación unitaria de compresión en el hormigón.
% % Si lo supera, entonces no se necesita elemento de borde.
% colorPalette = colormap(winter(2));
% colorPalette = colorPalette / max(colorPalette(:));
% figure2 = figure('InvertHardcopy','off','Color',[1 1 1],'position',[680 341 968 637]);
% axes2 = axes('Parent',figure2);
% hold on
% color_alea = colorPalette(2,:);
% plot([0; curvature(:,end)*100], [0; M(:,end)], 'linewidth', 3, 'Color', color_alea)
% plot([0; curvature_neg(:,end)*100], [0; M_neg(:,end)], 'linewidth', 3, 'Color', color_alea)
% xline(phi_u,'--r', '$$\phi_u$$', 'Interpreter', 'latex', 'fontsize', 20,'linewidth', 3)
% xline(-phi_u,'--r', '$$-\phi_u$$', 'Interpreter', 'latex', 'fontsize', 20, 'linewidth', 3)
% yline(0)
% xline(0)
% hold off
% xlabel('Curvature (phi) [1/m]')
% ylabel('Moment (M) [tonf]')
% grid on
% box on
% set(axes2, 'FontSize', 20);
% legend(sprintf('P_u = %.1f [tonf]', max(Pu_)));
% title('ec_{max} = 0.003')
% 
% % Hay un c_lim que me limita el máximo 'c' que puedo obtener, si este 'c'
% % supera c_lim, entonces se necesita confinar un mínimo largo de c_c = c-c_lim
% % (además hay otras consideraciones) que pueden agrandar este c_c
% fprintf('Elemento de Borde-------------------------------\n')
% c_lim = l_w/(600*(max([du/h_w; 0.007]))); % cm
% if c_pos >= c_lim 
%     fprintf('c_pos = %.1f [cm] > c_lim = %.1f[cm], hay que confinar arriba c_c = c-c_lim = %.0f [cm]\n',c_pos,c_lim, c_pos - c_lim)
% else
%     fprintf('c_pos = %.1f [cm] < c_lim = %.1f[cm], no hay que confinar arriba c_c = c-c_lim = %.0f [cm]\n',c_pos,c_lim, c_pos - c_lim)
% end
% 
% if c_neg_ >= c_lim 
%     fprintf('c_neg = %.1f [cm] > c_lim = %.1f[cm], hay que confinar abajo c_c = c-c_lim = %.0f [cm]\n',c_neg_, c_lim, c_neg_ - c_lim)
% else
%     fprintf('c_neg = %.1f [cm] < c_lim = %.1f[cm], no hay que confinar abajo c_c = c-c_lim = %.0f [cm]\n',c_neg_, c_lim, c_neg_ - c_lim)
% end
% 
% %% Diseño de elemento de borde
% 
% % altura del elemento de borde
% muvu = max(abs(Mu_./abs(internalForcesData(:,3)))); % metros
% h_be = max(l_w,muvu/4);


