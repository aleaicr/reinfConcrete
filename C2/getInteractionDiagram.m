function [Mn, Pn, phiMn, phiPn] = getInteractionDiagram(Section)
% Tarea 2 - Hormigón Armado Avanzado
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras R. - Gabriel Ramos V.
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
% es range
% es_vect = (es_min:(es_max-es_min)/(n_es-1):es_max).';
es_vect = [linspace(es_min, 0, 0.5*n_es).'; linspace(0, es_max, 0.5*n_es).'];
es_vect(es_vect==0) = []; % el código no calcula con este valor (probar nuevamente)
es_length = length(es_vect);

% negative section
Section_neg = Section;
Section_neg.h = flip(h);
Section_neg.b = flip(b);
Section_neg.as = flip(as); % cm2
Section_neg.d = sum(h) - flip(d); % cm
Section_neg.PC = sum(h) - PC; % cm

%% Interaction diagram computation
% init variables
Mn = zeros(es_length,1); % kgf
Pn = zeros(es_length,1); % kgf
phi_min = zeros(es_length,1);

Mn_neg = zeros(es_length,1);
Pn_neg = zeros(es_length,1);
phi_min_neg = zeros(es_length,1);

% Computation of values
for i = 1:es_length
    Section.es_val = es_vect(i);
    [Mn(i), Pn(i), phi_min(i)] = getMn_esBased(Section); % return en kgf y cm
    Section_neg.es_val = es_vect(i);
    [Mn_neg(i), Pn_neg(i), phi_min_neg(i)] = getMn_esBased(Section_neg); % return en kgf y cm
end

%% kgf, cm to tonf, m
Mn = Mn/1000/100; % tonf-m
Pn = Pn/1000; % tonf
phiMn = phi_min.*Mn;
phiPn = phi_min.*Pn;
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
plot(Mu_, Pu_, '.', 'color', '#C1330C')
grid on
xlabel('M_n [tonf-m]')
ylabel('P_n [tonf]')
box on
set(axe1,'FontSize',20);
hold off
annotation(figure1,'textbox',...
    [0.600297012302279 0.567185289957569 0.115871704745167 0.0650636492220651],...
    'Color',[0.466666666666667 0.674509803921569 0.188235294117647],...
    'String','\phi P_n, \phi M_n',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);
annotation(figure1,'textbox',...
    [0.715411247803157 0.780763790664782 0.115871704745167 0.0650636492220651],...
    'Color',[0.494117647058824 0.184313725490196 0.556862745098039],...
    'String','P_n, M_n',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);
annotation(figure1,'textbox',...
    [0.549330404217918 0.342291371994344 0.115871704745167 0.0650636492220651],...
    'Color',[0.635294117647059 0.0784313725490196 0.184313725490196],...
    'String','P_u, M_u',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);



end
