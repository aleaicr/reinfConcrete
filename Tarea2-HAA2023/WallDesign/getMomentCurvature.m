function [M, curvature, M_neg, curvature_neg] = getMomentCurvature(Section)
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
h = Section.h;
% r = Section.r;
% nBars = Section.nBars;
% diams = Section.diams;
% ecu = Section.ecu;
% nLayers = Section.nLayers;
% layers = Section.layers;
d = Section.d;
as = Section.as;
% P0 = Section.P0; % kgf
PC = Section.PC;
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
% define vector of ec range
ec_vect = (ec_min:(ec_max-ec_min)/(n_ec-1):ec_max).';
ec_length = length(ec_vect);

% init variables
M = zeros(ec_length,1);
M_neg = zeros(ec_length,1);
curvature = zeros(ec_length,1);
curvature_neg = zeros(ec_length,1);

% define negative section (flip the section)
Section_neg = Section;
Section_neg.as = flip(as); % cm2
Section_neg.d = sum(h) - flip(d); % cm
Section_neg.PC = sum(h) - PC; % cm

% compute moment and curvatures for all ec range
for i = 1:ec_length
    ec_val = ec_vect(i);
    [M(i), curvature(i)] = getMn_ecBased(Section, N, ec_val); % kgf, cm
    [M_neg(i),curvature_neg(i)] = getMn_ecBased(Section_neg,N,ec_val); % kgf, cm
end

% transform units and nagatives
M = M/1000/100; % tonf-m
M_neg = M_neg/1000/100; % tonf-m
curvature_neg = -curvature_neg; % 1/cm
M_neg = -M_neg; % tonf-m

%% Plot
figure1 = figure('InvertHardcopy','off','Color',[1 1 1],'position',[680   341   968   637]);
axes1 = axes('Parent',figure1);
plot([0; curvature], [0; M],'color','r','linewidth',3)
hold on
plot([0; curvature_neg], [0; M_neg],'color','r','linewidth',3)
% scatter(curvature(end),M(end),'^r','linewidth',10)
hold off
xlabel('Curvature (phi) [1/cm]')
ylabel('Momento (M) [tonf]')
grid on
box on
set(axes1,'FontSize',20);

end

