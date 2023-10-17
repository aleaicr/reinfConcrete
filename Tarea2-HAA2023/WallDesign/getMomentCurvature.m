function [M, curvature] = getMomentCurvature(Section)
% This function computes and plot the moment-curvature diagram by varying
% the 'ec' values and using the lowest Pu_ value
%
% Inputs
% Section: Section Properties
% Outputs
%
% Notes
%

%% Disarrange
% fc = Section.fc;
% fy = Section.fy;
% Es = Section.Es;
% b = Section.b;
% h = Section.h;
% r = Section.r;
% nBars = Section.nBars;
% diams = Section.diams;
% ecu = Section.ecu;
% nLayers = Section.nLayers;
% layers = Section.layers;
% d = Section.d;
% as = Section.as;
% P0 = Section.P0; % kgf
% PC = Section.PC;
% beta1_val = Section.beta1_val;
% es_min = Section.ess(1);
% es_max = Section.ess(2);
% n_es = Section.ess(3);
% Mu_ = Section.Mu_;
Pu_ = Section.Pu_;
ec_min = Section.ecc(1);
ec_max = Section.ecc(2);
n_ec = Section.ecc(3);

%% Previous
N = max(Pu_)*1000; %kgf

%% Get curvature
ec_vect = (ec_min:(ec_max-ec_min)/(n_ec-1):ec_max).';
ec_length = length(ec_vect);
M = zeros(ec_length,1);
curvature = zeros(ec_length,1);
for i = 1:ec_length
    ec_val = ec_vect(i);
    [M(i), curvature(i)] = getMn_ecBased(Section, N, ec_val); %kgf, cm
end

M = M/1000/100; % tonf-m

%% Plot
figure
plot(curvature, M)
xlabel('Curvature (phi) [1/m]')
ylabel('Momento (M) [tonf]')
grid on


end

