function [Mn, Pn, phiMn, phiPn, phi_curvature] = getInteractionDiagram(Section)
% Tarea 2 - Hormigón Armado Avanzado
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras R. - Gabriel Ramos V.
%
%%
% This function generates the Interaction Diagram (DI) and the 
% Moment-Curvature (M-C) diagrams of a reinforced concrete section
%
% INPUTS
% Section: (Struct) All the sections properties that are disarranged in the
% first section of this script
%
% OUTPUTS
% Mn: (double vector)
% Pn
% phiMn
% phiPn
% phi_curvature
%
% Notes
% 

%% Disarrange
fc = Section.fc;
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
P0 = Section.P0; % kgf
PC = Section.PC;
% beta1_val = Section.beta1_val;
es_min = Section.ess(1);
es_max = Section.ess(2);
n_es = Section.ess(3);
Mu_ = Section.Mu_;
Pu_ = Section.Pu_;

%% previous
es_vect = (es_min:(es_max-es_min)/(n_es-1):es_max).';
es_vect(es_vect==0) = [];
es_length = length(es_vect);
Section_neg = Section;
Section_neg.as = flip(as); % cm2
Section_neg.d = sum(h) - flip(d); % cm
Section_neg.PC = sum(h) - PC; % cm

%% Interaction diagram computation
% init variables
Mn = zeros(es_length,1); % kgf
Pn = zeros(es_length,1); % kgf
phi_min = zeros(es_length,1);
phi_curvature = zeros(es_length,1);

Mn_neg = zeros(es_length,1);
Pn_neg = zeros(es_length,1);
phi_min_neg = zeros(es_length,1);
phi_curvature_neg = zeros(es_length,1);

% positive side
for i = 1:length(es_vect)
    Section.es_val = es_vect(i);
    [Mn(i), Pn(i), phi_min(i), phi_curvature(i)] = getMn_esBased(Section); % return en kgf y cm
end

Mn = Mn/1000/100; % tonf-m
Pn = Pn/1000; % tonf
phiMn = phi_min.*Mn;
phiPn = phi_min.*Pn;

% Negative side
for i = 1:es_length
    Section_neg.es_val = es_vect(i);
    [Mn_neg(i), Pn_neg(i), phi_min_neg(i), phi_curvature_neg(i)] = getMn_esBased(Section_neg); % return en kgf y cm
end

Mn_neg = -Mn_neg/1000/100; % tonf-m
Pn_neg = Pn_neg/1000; % tonf
phiMn_neg = phi_min_neg.*Mn_neg; % tonf-m
phiPn_neg = phi_min_neg.*Pn_neg; % tonf

%% maximum phiPn correction
phiPn_max = 0.8*0.65*P0/1000; % tonf
phiPn = min(phiPn, phiPn_max); % tonf
phiPn_neg = min(phiPn_neg, phiPn_max); % tonf

%% Plot
% Interaction Diagram
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1]);
set(figure1,'Position',[496 271 1138 707])
axe1 = axes('Parent',figure1);
plot(Mn, Pn, 'color', '#460CB2', 'linewidth', 2)
hold on
yline(0.35*fc*sum(b.*h)/1000,'--r',"$$0.35f'_cA_g$$",'Interpreter', 'latex','fontsize',20)
plot(Mn_neg,Pn_neg,'color', '#460CB2', 'linewidth', 2)
plot(phiMn, phiPn,'color','#57A413', 'linewidth', 2)
plot(phiMn_neg, phiPn_neg,'color','#57A413', 'linewidth', 2)
plot(Mu_, Pu_, 'o', 'color', '#C1330C')
% plot(Mu2_, Pu_, 'o', 'color', '#C1330C')
grid on
xlabel('M_n [tonf-m]')
ylabel('P_n [tonf]')
box on
set(axe1,'FontSize',20);
hold off

% % Moment - Curvature
% figure
% plot(phi_curvature, Mn)
% hold on
% plot(phi_curvature_neg, Mn_neg)
% hold off
% xlabel('\phi curvature')
% ylabel('Mn [tonf]')

end
