%% Verificación Dimensionamiento FUNDACIONES
% Alexis Contreras R - Gabriel Ramos V.
% Universidad Técnica Federico Santa María - Campus SJ
% Hormigón Armado Avanzado - TAREA 3
%

%% Init
clear variables
close all
clc

%% Inputs
loadsfileName = 'loads.xlsx';
SAP2000ResultsFileName = 'ModeloconVF.xlsx';
H_total = 200;  % cm
B = 200; % cm
Largo = 200; % cm
gamma_ho = 2.5; % [tonf/m3]
despl_max = 0.4; % [mm]
sigma_sism_adm = 13; % [kgf/cm2]
a = 0.25;

COLUMNSIDS = [1; 2; 3; 4; 5; 6; 7; 8; 9];

%% Read load data
loads = readtable(loadsfileName);

% Filtrar Datos
D = loads(strcmp(loads.OutputCase, 'D'), :);
L = loads(strcmp(loads.OutputCase, 'L'), :);
SO = loads(strcmp(loads.OutputCase, 'SO'), :);
Ex = loads(strcmp(loads.OutputCase, 'Ex'), :);
Ey = loads(strcmp(loads.OutputCase, 'Ey'), :);
Ez = loads(strcmp(loads.OutputCase, 'Ez'), :);
D = D{:,3:7};
L = L{:,3:7};
SO = SO{:,3:7};
Ex = Ex{:,3:7};
Ey = Ey{:,3:7};
Ez = Ez{:,3:7};

% Name of the variables
variableNames = loads.Properties.VariableNames;

%% Load combinations
% Simultaneidad sísmica
signs = [-1, 1]; 
E1 = zeros(size(D,1), size(D,2),6);
E2 = E1;
E3 = E1;
E4 = E1;
l = 0;
for i = 1:length(signs)
    for j = 1:length(signs)
        for k = 1:length(signs)
            l = l + 1;
            % Calcular E1, E2, E3, E4 con combinaciones de signos
            E1(:,:,l) = signs(i)*Ex + signs(j)*0.3*Ey + signs(k)*0.6*Ez;
            E2(:,:,l) = signs(i)*0.3*Ex + signs(j)*Ey + signs(k)*0.6*Ez;
            E3(:,:,l) = signs(i)*0.6*Ex + signs(j)*0.2*Ey + signs(k)*Ez;
            E4(:,:,l) = signs(i)*0.2*Ex + signs(j)*0.6*Ey + signs(k)*Ez;
        end
    end
end
E_all = cat(3, E1, E2, E3, E4); % concatenar

% Combinaciones de carga
ASD1 = zeros(size(D,1), size(D,2),l*4);
ASD2 = ASD1;
LRFD1 = ASD1;
LRFD2 = ASD1;
for i = 1:(l*4)
    E = E_all(:,:,l);
    ASD1(:,:,i) = D + 0.75*a*L + 0.75*SO + E;
    ASD2(:,:,i) = D + E;
    LRFD1(:,:,i) = 1.2*D + a*L + SO + 1.4*E;
    LRFD2(:,:,i) = 0.9*D + 1.4*E; 
end
ASD_all = cat(3, ASD1, ASD2);
LRFD_all = cat(3, LRFD1, LRFD2);

%% Dimensionamiento
% SOIL PRESSURE
% Import Soil Pressure Data - UNITS: TONF M 
SheetToImport = "Element Soil Pres - Area Shells";                          % Sheet that have the soil data
soilDataTable = readtable(SAP2000ResultsFileName, 'Sheet', SheetToImport);  % Import soil data from SAP2000 results excel file

% Maximum soil preassure for ASD combinations
ASDcombsTable = soilDataTable(contains(soilDataTable.OutputCase, 'ASD'), :); % Table with only ASD combinations (must be defined with ASD in the load case name)
pres_max = -min(ASDcombsTable.Pressure);                                    % Maximum soil pressure
ASD_combo_max_pres = ASDcombsTable(ASDcombsTable.Pressure == -pres_max, :); % 

% print
fprintf('Presión máxima del suelo--------\n')
if pres_max/10 < sigma_sism_adm
    fprintf('P_max = %.2f[kgf/cm2] < sigma_sism_adm = %.2f[kgf/cm2] OK\n', pres_max/10, sigma_sism_adm)
else
    fprintf('P_max = %.2f[kgf/cm2] > sigma_sism_adm = %.2f[kgf/cm2] NO OK\n', pres_max/10, sigma_sism_adm)
end
fprintf('Ocurre para la combinación %s\n', string(ASD_combo_max_pres.OutputCase(1)))

% MAXIMUM VERTICAL JOINT DISPLACEMENT
% Import vertical displacement data - UNITS: M
SheetToImport = "Joint Displacements";
jointTable = readtable(SAP2000ResultsFileName, 'Sheet', SheetToImport);

% Maximum joint displacement (in the entire foundation)
ASDcombsTable = jointTable(contains(jointTable.OutputCase, 'ASD'), :);
[max_disp, index] = max(abs(ASDcombsTable.U3));
max_disp_sign = sign(ASDcombsTable.U3(index)); % m
ASD_combo_max_disp = ASDcombsTable(abs(ASDcombsTable.U3) == max_disp, :);

% print
fprintf('\nDesplazamiento máximo del suelo--------\n')
if max_disp*1000 <= despl_max
    fprintf('Despl_max = %.3f[mm] < %.3f[mm] OK\n', max_disp_sign*max_disp*1000, despl_max)
else
    fprintf('Despl_max = %.3f[mm] > %.3f[mm] NO OK\n', max_disp_sign*max_disp*1000, despl_max)
end
fprintf('Ocurre para la combinación %s\n',string(ASD_combo_max_disp.OutputCase(1)))


% MOMENTO VOLCANTE
% Import 


% DESLIZAMIENTO LATERAL




% % Inicializando variables
% Mx = zeros(10,64);
% N = Mx;
% e = Mx;
% My = Mx;
% 
% Mx = ASD_all(:,4) + ASD_all(:,2)*H;
% My = ASD_all(:,5) + ASD_all(:,3)*H;
% WF = 
% N = 

%% Design
% GET MAXIMUM M11+, M11-, M22+, M22- IN SHELL ELEMENTS - TONF-M
% Import vertical displacement data - UNITS: M
SheetToImport = "Element Forces - Area Shells";
shellTable = readtable(SAP2000ResultsFileName, 'Sheet', SheetToImport);

% Maximum joint displacement (in the entire foundation)
LRFDcombsTable = shellTable(contains(shellTable.OutputCase, 'LRFD'), :);
LRFDcombs = LRFDcombsTable.OutputCase;
LRFDcombs_names = unique(LRFDcombs);
max_m11combs = zeros(length(LRFDcombs_names),1);
min_m11combs = max_m11combs;
max_m22combs = max_m11combs;
min_m22combs = max_m11combs;
fprintf('Loop in the table (for each LRFD comb)\n')
for i = 1:length(LRFDcombs_names)
    % Create a table for this combination only
    LRFDTable_aux = LRFDcombsTable(contains(LRFDcombsTable.OutputCase, char(LRFDcombs_names(i))), :);
    
    % Save max MXX+- of each combination
    max_m11combs(i,1) = max(LRFDTable_aux.M11);
    min_m11combs(i,1) = min(LRFDTable_aux.M11);
    max_m22combs(i,1) = max(LRFDTable_aux.M22);
    min_m22combs(i,1) = min(LRFDTable_aux.M22);
end

% Get maximum MX+-
[max_m11, index_max11] = max(max_m11combs);
[min_m11, index_min11] = min(min_m11combs);
[max_m22, index_max22] = max(max_m22combs);
[min_m22, index_min22] = min(min_m22combs);

% Get combos where MX+- occurs
max11_combo_name = char(LRFDcombs_names(index_max11));
min11_combo_name = char(LRFDcombs_names(index_min11));
max22_combo_name = char(LRFDcombs_names(index_max22));
min22_combo_name = char(LRFDcombs_names(index_min22));

fprintf('Solicitaciones FLEXIÓN LRFD máximas fundación---------\n')
fprintf('maxM11 = %.2f[tonf-m]\n', max_m11)
fprintf('Ocurrió para la comb %s\n\n', max11_combo_name)
fprintf('minM11 = %.2f[tonf-m]\n', min_m11)
fprintf('Ocurrió para la comb %s\n\n', min11_combo_name)
fprintf('maxM22 = %.2f[tonf-m]\n', max_m22)
fprintf('Ocurrió para la comb %s\n\n', max22_combo_name)
fprintf('minM22 = %.2f[tonf-m]\n', min_m22)
fprintf('Ocurrió para la comb %s\n\n', min22_combo_name)


% GET M, N, V IN COLUMN
SheetToImport = "Element Forces - Frames";
% -> load data

% -> filter by column IDS (are not the only frame elements)

% -> for each column ID get the internal forces matrix (N,V2,V3,M2,M3)

% -> armar con DI



