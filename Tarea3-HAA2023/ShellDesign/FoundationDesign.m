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

%% Foundation Inputs
% SAP2000 results file Name
SAP2000fileDir = '..\ModeloconVF.xlsx'; 

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

% Reinforcement
s = 10; % cm                                                                % Spacing
diams = 12;  % mm                                                           % Reinforcement diameter
nBars = (b - s)/s + 1;                                                      % Number of bars
d = h - r; % cm                                                             % Effective depth

% ecu
ecu = 0.003;

% Minimum as/ag
rho_min = 1.8/1000;

% Bottom columns nodes (shear is too high there)
cNodes = [18;16;17;26;19;25;20;14;15;24];                                   % Vector with all the nodes that directly connect a frame and shear.

%% Previous Calculations
% Neglecting the axial load in the flexural design of the shell element
loadComb = 'LRFD';
[M11, M22, V13, V23] = getShellLoads(SAP2000fileDir, loadComb, cNodes);
minM11 = M11.min;
maxM11 = M11.max;
minM22 = M22.min;
maxM22 = M22.max;
maxV13 = V13.max;
maxV23 = V23.max;

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
% Section.Pu_ = 0;
% Section.n_N = n_N;
% Section.part = part;

%% FLEXURAL DESIGN
%% Print internal loads
fprintf('M11pos = %.2f en %s\n', maxM11, M11.max_OutputCase{1})
fprintf('M11neg = %.2f en %s\n', minM11, M11.min_OutputCase{1})
fprintf('M22pos = %.2f en %s\n', maxM22, M22.max_OutputCase{1})
fprintf('M22neg = %.2f en %s\n\n', minM22, M22.min_OutputCase{1})

%% M11
% get Mnpos11
Section.N = 0; % kgf
Section.Mu_ = maxM11;
[Mn11pos, phiMn11pos, ~, ~, ~, ~, ~, ~] = getMn(Section); % kgf - cm (inputs y otuputs)
Mn11pos = Mn11pos/1000/100; % tonf-m/m                                              % Porque b = 100cm = 1m
phiMn11pos = phiMn11pos/1000/100; % tonf-m/m

fprintf('phi%.0f@%.0f\n',diams, s)
fprintf('Mn_pos_11 = %.3f [tonf-m/m]\n', Mn11pos)
fprintf('phiMn_pos_11 = %.3f [tonf-m/m]\n', phiMn11pos)

% get Mnneg11
s = 10; % cm                                                                % Spacing
diams = 12;  % mm                                                          % Reinforcement diameter
nBars = (b - s)/s + 1;                                                      % Number of bars
as_types = 0.25*pi*(nBars.*(diams.'/10).^2); % cm^2
as = sum(as_types,2); % cm^2                                                % Vector of as in each layer of reinforcement
P0 = 0.85*fc*(ag-sum(as)) + sum(as*fy); % kgf
PC = (0.85*fc*sum(b.*h.*h_centroids) + sum(as.*d*(fy-0.85*fc)))/P0; % cm

Section2 = Section;
Section2.as = as;
Section2.P0 = P0;
Section2.PC = PC;

[Mn11neg, phiMn11neg, ~, ~, ~, ~, ~, ~] = getMn(Section2); % kgf - cm (inputs y otuputs)
Mn11neg = Mn11neg/1000/100; % tonf-m/m                                              % Porque b = 100cm = 1m
phiMn11neg = phiMn11neg/1000/100; % tonf-m/m

fprintf('phi%.0f@%.0f\n',diams, s)
fprintf('Mn_neg_11 = %.3f [tonf-m/m]\n', Mn11neg)
fprintf('phiMn_neg_11 = %.3f [tonf-m/m]\n', phiMn11neg)
% Cuantía
rho_disp = sum(Section.as + Section2.as)/ag;
if rho_disp>rho_min
    fprintf('rho_disp = %.4f > rho_min = %.4f OK\n\n', rho_disp, rho_min)
else
    fprintf('rho_disp = %.4f < rho_min = %.4f NO OK\n\n', rho_disp, rho_min)
end

%% M22
% get Mnpos22
s = 10; % cm                                                                % Spacing
diams = 16;  % mm                                                           % Reinforcement diameter
nBars = (b - s)/s + 1;                                                      % Number of bars
d = h - r; % cm                                                             % Effective depth
as_types = 0.25*pi*(nBars.*(diams.'/10).^2); % cm^2
as = sum(as_types,2); % cm^2                                                % Vector of as in each layer of reinforcement
P0 = 0.85*fc*(ag-sum(as)) + sum(as*fy); % kgf
aux = h;
aux(end) = [];
aux = [0; cumsum(aux)];
h_centroids = h/2 + aux; % cm
PC = (0.85*fc*sum(b.*h.*h_centroids) + sum(as.*d*(fy-0.85*fc)))/P0; % cm

Section3 = Section;
Section3.s = s;
Section3.diams = diams;
Section3.nBars = nBars;
Section3.d = d;
Section3.Mu_ = maxM22;
Section3.as = as;
Section3.PC = PC;

[Mn22pos, phiMn22pos, ~, ~, ~, ~, ~, ~] = getMn(Section3); % kgf - cm (inputs y otuputs)
Mn22pos = Mn22pos/1000/100; % tonf-m/m                                              % Porque b = 100cm = 1m
phiMn22pos = phiMn22pos/1000/100; % tonf-m/m

fprintf('phi%.0f@%.0f\n',diams, s)
fprintf('Mn_pos_22 = %.3f [tonf-m/m]\n', Mn22pos)
fprintf('phiMn_pos_22 = %.3f [tonf-m/m]\n', phiMn22pos)

% get Mnneg22
s = 10; % cm                                                                % Spacing
diams = 12;  % mm                                                           % Reinforcement diameter
nBars = (b - s)/s + 1;                                                      % Number of bars
as_types = 0.25*pi*(nBars.*(diams.'/10).^2); % cm^2
as = sum(as_types,2); % cm^2                                                % Vector of as in each layer of reinforcement
P0 = 0.85*fc*(ag-sum(as)) + sum(as*fy); % kgf
PC = (0.85*fc*sum(b.*h.*h_centroids) + sum(as.*d*(fy-0.85*fc)))/P0; % cm

Section4 = Section;
Section4.s = s;
Section4.diams = diams;
Section4.nBars = nBars;
Section4.d = d;
Section4.Mu_ = minM22;
Section4.as = as;
Section4.PC = PC;

[Mn22neg, phiMn22neg, ~, ~, ~, ~, ~, ~] = getMn(Section4); % kgf - cm (inputs y otuputs)
Mn22neg = Mn22neg/1000/100; % tonf-m/m                                              % Porque b = 100cm = 1m
phiMn22neg = phiMn22neg/1000/100; % tonf-m/m

fprintf('phi%.0f@%.0f\n',diams, s)
fprintf('Mn_pos_22 = %.3f [tonf-m/m]\n', Mn22neg)
fprintf('phiMn_pos_22 = %.3f [tonf-m/m]\n', phiMn22neg)
% Cuantía
rho_disp = sum(Section3.as + Section4.as)/ag;
if rho_disp>rho_min
    fprintf('rho_disp = %.4f > rho_min = %.4f OK\n\n', rho_disp, rho_min)
else
    fprintf('rho_disp = %.4f < rho_min = %.4f NO OK\n\n', rho_disp, rho_min)
end

%% SHEAR DESIGN
% Print internal loads
fprintf('maxV13 = %.2f en %s\n', maxV13, V13.max_OutputCase{1})
fprintf('maxV23 = %.2f en %s\n\n', maxV23, V23.max_OutputCase{1})

% Check concrete's section shear capacity ONE-DIRECTION
fprintf('ONE-DIRECTION\n')
Vc = 0.53*sqrt(fc)*b*(h-r)/1000; % tonf (per meter)
phiMin_corte = 0.75;
Vu_Corteviga = 40.15;      % Manualmente
if phiMin_corte*Vc > Vu_Corteviga
    fprintf('phiVc = %.2f [tonf] (por metro) > maxVu = %.2f OK\n\n', phiMin_corte*Vc, Vu_Corteviga)
else
    fprintf('phiVc = %.2f [tonf] (por metro) < maxVu = %.2f NO OK\n\n', phiMin_corte*Vc, Vu_Corteviga)
end
if phiMin_corte*Vc > Vu_Corteviga
    fprintf('phiVc = %.2f [tonf] (por metro) > maxVu = %.2f OK\n\n', phiMin_corte*Vc, Vu_Corteviga)
else
    fprintf('phiVc = %.2f [tonf] (por metro) < maxVu = %.2f NO OK\n\n', phiMin_corte*Vc, Vu_Corteviga)
end

% Check punching shear capacity TWO-DIRECTIONS
fprintf('TWO-DIRECTIONS\n')
% Beams 65x65
Vu = 84.16;
b0 = 123*4; % cm
beta = 1; % 65/65
alfa_s = 20; % Worst case (border corner)
d = (65-7); % cm
Vc1 = 1.0*sqrt(fc)*b0*(d); % kgf
Vc2 = 0.53*(1+2/beta)*sqrt(fc)*b0*d; % kgf
Vc3 = 0.27*(alfa_s*d/b0 + 2)*sqrt(fc)*b0*d; % kgf
Vc_65x65 = min([Vc1; Vc2; Vc3])/1000; % tonf

if phiMin_corte*Vc_65x65 > Vu
    fprintf('phiVc_65x65 = %.2f[tonf] > maxVu = %.2f [tonf] OK \n\n', phiMin_corte*Vc_65x65, Vu)
else
    fprintf('phiVc_65x65 = %.2f[tonf] < maxVu = %.2f [tonf] NO OK\n\n', phiMin_corte*Vc_65x65, Vu)
end
% Beams 70x70
Vu = 72.28;
b0 = 133*4; % m
beta = 1; % 70/70
alfa_s = 30; % Worst case (border column)
d = 70-7;
Vc1 = 1.0*sqrt(fc)*b0*d; % kgf
Vc2 = 0.53*(1+2/beta)*sqrt(fc)*b0*d; % kgf
Vc3 = 0.27*(alfa_s*d/b0 + 2)*sqrt(fc)*b0*d; % kgf
Vc_70x70 =  min([Vc1; Vc2; Vc3])/1000; % tonf 
if phiMin_corte*Vc_70x70 > Vu
    fprintf('phiVc_70x70 = %.2f[tonf] > maxVu = %.2f [tonf] OK\n\n', phiMin_corte*Vc_70x70, Vu)
else
    fprintf('phiVc_70x70 = %.2f[tonf] < maxVu = %.2f [tonf] NO OK\n\n', phiMin_corte*Vc_70x70, Vu)
end


