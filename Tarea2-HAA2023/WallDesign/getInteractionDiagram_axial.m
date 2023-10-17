function [Mn, Pn, phiMn, phiPn] = getInteractionDiagram_axial(Section)
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
fc = Section.fc; % kgf/cm2
fy = Section.fy; % kgf/cm2
% Es = Section.Es;
b = Section.b; % cm
h = Section.h; % cm
% r = Section.r;
% nBars = Section.nBars;
% diams = Section.diams;
% ecu = Section.ecu;
% nLayers = Section.nLayers;
% layers = Section.layers;
d = Section.d; % cm
as = Section.as; %cm2
P0 = Section.P0; % kgf
PC = Section.PC; % cm
% beta1_val = Section.beta1_val;
% es_min = Section.ess(1);
% es_max = Section.ess(2);
% n_es = Section.ess(3);
Mu_ = Section.Mu_;
Pu_ = Section.Pu_;
n_N = Section.n_N;

%% previous
Section_neg = Section;
Section_neg.as = flip(as); % cm2
Section_neg.d = sum(h) - flip(d); % cm
Section_neg.PC = sum(h) - PC; % cm

% compresión pura
Pcompr = P0; %kgf

% tracción pura
Ptracc = -sum(as)*fy; % kgf

% axial load range
% N_step = (Pcompr-Ptracc)/(n_N-1);
N_vect = (Ptracc:(Pcompr-Ptracc)/(n_N-1):Pcompr).'; %kgf
N_vect(N_vect==0) = [];
N_length = length(N_vect); % debería ser n_N + 2
Pn = N_vect; % kgf
Pn_neg = N_vect; %kgf

%% Interaction diagram computation
% init variables
Mn = zeros(N_length,1); % kgf-cm
phi_min = zeros(N_length,1);
Mn_neg = zeros(N_length,1); %kgf-cm
phi_min_neg = zeros(N_length,1);

% computation of Mn in function of axial loading
for i = 1:N_length
    Section.N = N_vect(i); % kgf
    Section_neg.N = N_vect(i); % kgf
    disp(N_vect(i))
    [Mn(i), phi_min(i)] = getMn_axial(Section); % return en kgf y cm
    disp(N_vect(i))
    [Mn_neg(i), phi_min_neg(i)] = getMn_axial(Section_neg); % return en kgf y cm
end

% Cambio de unidades
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
yline(phiPn_max,'--r',"$$\phi P_{n,max}$$",'Interpreter', 'latex')
yline(0.35*fc*sum(b.*h)/1000,'--r',"$$0.35f'_cA_g$$",'Interpreter', 'latex')
plot(Mn_neg, Pn_neg,'color', '#460CB2', 'linewidth', 2)
plot(phiMn, phiPn,'color','#57A413', 'linewidth', 2)
plot(phiMn_neg, phiPn_neg,'color','#57A413', 'linewidth', 2)
plot(Mu_, Pu_, 'o', 'color', '#C1330C')
grid on
xlabel('M_n [tonf-m]')
ylabel('P_n [tonf]')
% title('Interaction Diagram')
% legend('R. Nominal','','R .Diseño','','Solicitaciones')
% Set the remaining axes properties
annotation('textbox',...
    [0.5638078817734 0.392857142857151 0.107374384236453 0.0761904761904769],...
    'Color',[0.635294117647059 0.0784313725490196 0.184313725490196],...
    'String','$$M_u, P_u$$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off');
annotation('textbox',...
    [0.704625041050904 0.489717261904767 0.107374384236453 0.076190476190477],...
    'Color',[0.466666666666667 0.674509803921569 0.188235294117647],...
    'String','$$\phi M_n, \phi P_n$$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off');
annotation('textbox',...
    [0.757100123152709 0.668556547619051 0.107374384236453 0.0761904761904766],...
    'Color',[0.286274509803922 0 0.431372549019608],...
    'String','$$M_n, P_n$$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off');
box on
set(axe1,'FontSize',20);
hold off

end
