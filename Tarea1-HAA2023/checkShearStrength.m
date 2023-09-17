%% Shear

%% Init
clear variables
close all
clc

%% Inputs
% Loads
Vu = 50; % tonf
Vsism = 60; % tonf % Si el corte es mayor a 

% Meterials
fc = 300; % kgf/cm^2
fy = 4200; % kgf/cm^2
Es = 2.1*10^6; %kgf/cm^2

% Section
b = 40; % cm
h = 60; % cm
r = 5; % cm
d = 55; % cm

% Min. Factor
phi = 0.75; % Capacity Design, phi = 0.75, Struct.Analysis Design, phi = 0.6;

% Configuración
Estribos = 2; % Estribo simple 1, estribo doble 2
Trabas = 1; % Número de trabas

%% Concrete
if Vsims>Vu/2 && max(abs(Pu_))<= (b*h)*fc/20
    Vc = 0;
else
    Vc = 0.53*sqrt(fc)*b*d;
end

%% Steel
Vs_req = Vu/phi - Vc;
Av_s_req = Vs_req*1000/(fy*d);
Av_s_min = max(0.25*sqrt(fc)*b/fy,3.5*b/fy);

%% Configuración
n_barras = 2*Estribos + Trabas;

%% Av_s_Buscar
Av_s_buscar = max(Av_s_req*100/n_barras,av_s_min*100/n_barras);

