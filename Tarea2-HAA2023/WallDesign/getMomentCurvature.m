function [M, curvature, M_neg, curvature_neg, c, c_neg] = getMomentCurvature(Section, N_partitions)
% Tarea 2 - Hormigón Armado Avanzado
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras R. - Gabriel Ramos V.
%%
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
b = Section.b;
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
es_min_2 = Section.ess_2(1);
es_max_2 = Section.ess_2(2);
n_es_2 = Section.ess_2(3);
% Mu_ = Section.Mu_;
N = Section.Pu_;
ec_min = Section.ecc(1);
ec_max = Section.ecc(2);
n_ec = Section.ecc(3);

%% Previous
N_vect = linspace(min(N), max(N), N_partitions)*1000; % kgf

%% Get curvature
% define vector of ec range
ec_vect = logspace(log10(ec_min),log10(ec_max),n_ec);
% ec_vect = (ec_min:(ec_max-ec_min)/(n_ec-1):ec_max).';
ec_length = length(ec_vect);
es_vect = [linspace(es_min_2, 0, 0.5*n_es_2).'; linspace(0, es_max_2, 0.5*n_es_2).'];
es_vect(es_vect==0) = [];

% init variables
M = zeros(ec_length,N_partitions);
M_neg = zeros(ec_length,N_partitions);
curvature = zeros(ec_length,N_partitions);
curvature_neg = zeros(ec_length,N_partitions);
c = zeros(ec_length,N_partitions);
c_neg = zeros(ec_length,N_partitions);

% define negative section (flip the section)
Section_neg = Section;
Section_neg.h = flip(h);
Section_neg.b = flip(b);
Section_neg.as = flip(as); % cm2
Section_neg.d = sum(h) - flip(d); % cm
Section_neg.PC = sum(h) - PC; % cm

% compute moment and curvatures for all ec range
for i = 1:ec_length
    for j = 1:length(N_vect)
        ec = ec_vect(i);
        [M(i,j), curvature(i,j), c(i,j)] = getMn_ecBased(Section, N_vect(j), ec, es_vect); % kgf, cm
        [M_neg(i,j), curvature_neg(i,j), c_neg(i,j)] = getMn_ecBased(Section_neg, N_vect(j), ec, es_vect); % kgf, cm
        fprintf('%.0f, %.0f\n',i,j)
    end
end

% kgf, cm to tonf, m
M = M/1000/100; % tonf-m
M_neg = -M_neg/1000/100; % tonf-m
curvature_neg = -curvature_neg; % 1/cm

%% Plot
colorPalette = colormap(winter(N_partitions));
colorPalette = colorPalette / max(colorPalette(:));
figure1 = figure('InvertHardcopy','off','Color',[1 1 1],'position',[680 341 968 637]);
axes1 = axes('Parent',figure1);
hold on
legends = cell(2*N_partitions, 1);
for i = 1:N_partitions
    color_alea = colorPalette(i,:);
    plot([0; curvature(:,i)], [0; M(:,i)], 'linewidth', 3, 'Color', color_alea)
    plot([0; curvature_neg(:,i)], [0; M_neg(:,i)], 'linewidth', 3, 'Color', color_alea)
    legends{2*i-1} = ['N = ' num2str(N_vect(i)/1000) ' [tonf]'];
    legends{2*i} = '';
end
hold off
xlabel('Curvature (phi) [1/cm]')
ylabel('Moment (M) [tonf]')
grid on
box on
set(axes1, 'FontSize', 20);
legend(cellstr(legends));

end

