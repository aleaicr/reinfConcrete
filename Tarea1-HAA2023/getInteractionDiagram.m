function [Mn, Pn, phiMn, phiPn] = getInteractionDiagram(Section)
% This function generates the interaction diagram

fc = Section.fc;
fy = Section.fy;
Es = Section.Es;
b = Section.b;
h = Section.h;
% r = Section.r;
% nBars = Section.nBars;
% diams = Section.diams;
ecu = Section.ecu;
% nLayers = Section.nLayers;
% layers = Section.layers;
d = Section.d;
as = Section.as;
P0 = Section.P0; % kgf
PC = Section.PC;
beta1_val = Section.beta1_val;
es_min = Section.ess(1);
es_max = Section.ess(2);
n_es = Section.ess(3);
Mu_ = Section.Mu_;
Pu_ = Section.Pu_;

%% Compresión Pura
Pn_c_pura = P0/1000; % tonf
phiPn_c_Pura = Pn_c_pura*0.65; % tonf

%% Tracción Pura
Pn_t_pura = -sum(as)*fy/1000; % tonf
phiPn_t_pura = 0.9*Pn_t_pura; % tonf

%% Flexión Pura
Section.Pu = 0; % tonf
[Mn_pura, phiMn_pura, ~, ~, ~, ~, ~, ~]= getMn(Section); % return en kgf - cm
Mn_pura = Mn_pura/1000/100; % tonf
phiMn_pura = phiMn_pura/1000/100; %tonf

%% Falla Balanceada
es_b = 0.002;
d_last = max(d); % cm
c_b = d_last*ecu/(ecu+es_b); % cm
ess_b = (c_b-d)/c_b*ecu; 
phi_b = phi(min(ess_b));
sigma_steel_b = sigmaReinf(ess_b,fy,Es); % kgf/cm2
Fsteel_b = as.*sigma_steel_b; % kgf
Cc_b = 0.85*fc*beta1_val*c_b*b; % kgf
Mn_b = (Cc_b*(PC-beta1_val*c_b/2) + sum(Fsteel_b.*(PC-d)))/1000/100; % tonf
phiMn_b = Mn_b*phi_b; %tonf
Pn_b = (Cc_b + sum(Fsteel_b))/1000; % tonf
phiPn_b = phi_b*Pn_b; % tonf

%% For Loop
es_vect = (es_min:(es_max-es_min)/(n_es-1):es_max).';
Mn = zeros(length(es_vect),1);
Pn = Mn;
phi_min = Mn;
for i = 1:length(es_vect)
    Section.es_eval = es_vect(i);
    [Mn(i), Pn(i), phi_min(i)] = getMn_esBased(Section); % return en kgf y cm
end
Mn = Mn/1000/100; % tonf-m
Pn = Pn/1000; % tonf
phiMn = phi_min.*Mn;
phiPn = phi_min.*Pn;

%% Juntarlas
% CP,TP,FP,FB
markedMn = [0;0;Mn_pura;Mn_b]; % tonf-m
markedPn = [Pn_c_pura; Pn_t_pura; 0; Pn_b]; % tonf
markedphiMn = [0; 0; phiMn_pura; phiMn_b]; % tonf-m
markedphiPn = [phiPn_c_Pura; phiPn_t_pura; 0; phiPn_b]; %tonf

%% Corregir phiPn por phiPn_max
phiPn_max = 0.8*0.65*P0/1000; % tonf
phiPn = min(phiPn, phiPn_max); % tonf
markedphiPn = min(markedphiPn, phiPn_max); % tonf

%% Para el lado negativo
Section_neg = Section;
Section_neg.as = flip(as);
Section_neg.d = h - flip(d);
Section_neg.PC = h - PC;

es_vect = (es_min:(es_max-es_min)/(n_es-1):es_max).';
Mn_neg = zeros(length(es_vect),1);
Pn_neg = Mn;
phi_min_neg = Mn;
for i = 1:length(es_vect)
    Section_neg.es_eval = es_vect(i);
    [Mn_neg(i), Pn_neg(i), phi_min_neg(i)] = getMn_esBased(Section_neg); % return en kgf y cm
end
Mn_neg = - Mn_neg/1000/100; % tonf-m
Pn_neg = Pn_neg/1000; % tonf
phiMn_neg = phi_min_neg.*Mn_neg;
phiPn_neg = phi_min_neg.*Pn_neg;

phiPn_neg = min(phiPn_neg, phiPn_max); % tonf

%% Graficar
figure
plot(Mn, Pn, 'color', '#460CB2', 'linewidth', 2)
hold on
plot(phiMn, phiPn,'color','#57A413', 'linewidth', 2)
plot(Mu_, Pu_, 'o', 'color', '#C1330C')
hold off
grid on
xlabel('Mn')
ylabel('Pn')
title('Interaction Diagram')
legend('R. Nominal','R .Diseño','Solicitaciones')

% figure
% vect_Mn_b = (0:1:phiMn_b).';
% vect_Mn_b_length = length(vect_Mn_b);
% plot(Mn, Pn, 'color', '#460CB2', 'linewidth', 2)
% hold on
% plot(phiMn, phiPn,'color','#57A413', 'linewidth', 2)
% % plot(markedMn, markedPn,'o','color','r')
% % plot(markedphiMn, markedphiPn,'o','color','b')
% plot(Mu_, Pu_, 'o', 'color', '#C1330C')
% plot(vect_Mn_b, phiPn_b*ones(vect_Mn_b_length,1))
% hold off
% grid on
% xlabel('Mn')
% ylabel('Pn')
% title('Interaction Diagram')
% legend('R. Nominal','R .Diseño','Solicitaciones','Falla Balanceada')

figure
plot(Mn, Pn, 'color', '#460CB2', 'linewidth', 2)
hold on
plot(Mn_neg,Pn_neg,'color', '#460CB2', 'linewidth', 2)
plot(phiMn, phiPn,'color','#57A413', 'linewidth', 2)
plot(phiMn_neg, phiPn_neg,'color','#57A413', 'linewidth', 2)
plot(Mu_, Pu_, 'o', 'color', '#C1330C')
hold off
grid on
xlabel('Mn')
ylabel('Pn')
title('Interaction Diagram')
legend('R. Nominal + ','R. Nominal -','R .Diseño + ','R. Diseño -','Solicitaciones')

end
