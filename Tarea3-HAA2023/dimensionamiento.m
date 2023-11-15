

%% Init
clear variables
close all
clc

%% Inputs
fileName = 'datosTarea3.xlsx';
a = 0.25;
H_total = 200;  % del enunciado
B = 200; % 2m del dibujo
Largo = 200; % 2m este es el que tenemos que variar
gamma_ho = 2.5;

%% Read data
data = readtable(fileName);

% Filtrar Datos
D = data(strcmp(data.OutputCase, 'D'), :);
L = data(strcmp(data.OutputCase, 'L'), :);
SO = data(strcmp(data.OutputCase, 'SO'), :);
Ex = data(strcmp(data.OutputCase, 'Ex'), :);
Ey = data(strcmp(data.OutputCase, 'Ey'), :);
Ez = data(strcmp(data.OutputCase, 'Ez'), :);
D = D{:,3:7};
L = L{:,3:7};
SO = SO{:,3:7};
Ex = Ex{:,3:7};
Ey = Ey{:,3:7};
Ez = Ez{:,3:7};

% Nombre de las variables
variableNames = data.Properties.VariableNames;

%% Combinaciones de carga
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
for i = 1:(l*4)
    E = E_all(:,:,l);
    ASD1(:,:,i) = D + 0.75*a*L + 0.75*SO + E;
    ASD2(:,:,i) = D + E;
end
ASD_all = cat(3,ASD1,ASD2);


%% Dimensionamiento
% Inicializando variables
Mx = zeros(10,64);
N = Mx;
e = Mx;
My = Mx;

Mx = ASD_all(:,4) + ASD_all(:,2)*H;
My = ASD_all(:,5) + ASD_all(:,3)*H;
WF = 
N = 








