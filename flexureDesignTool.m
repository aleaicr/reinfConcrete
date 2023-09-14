%% Flexure Design Tool

%% Init
clear variables
close all
clc

%% Inputs
% Materials
fc = 300; % kgf/cm2
fy = 4200; % kgf/cm^2
Es = 2.1*10^6; %kgf/cm^2

% Section geometry
b = 30; % cm
h = 60; % cm
r = 5; % cm

% Reinforcement
nBars = [3; 3; 3];                                                          % Number of bars
diams = [2.2; 2.2; 2.2]; % cm

% Axial loading
Pu = 0; % tonf

% ecu
ecu = 0.003;

%% Previous Calculations
Ec = 15100*sqrt(fc); % kgf/cm^2                                                        % Moduluas of elasticity of unconfined concrete
nLayers = length(diams);                                                    % Number of Layers of longitudinal reinforcement
layers = (1:1:nLayers).';                                                   % Layer IDs
d = r + (h-2*r)/(nLayers-1)*(layers-1); % cm                                % Depth of the layers

%% Flexural Strength
% Area of steel
as = nBars*pi.*(diams/2).^2;    % cm^2                                      % vector (each row is a layer)

% Axial strength of the section
P0 = 0.85*fc*b*h+sum(as)*(fy - fc); % kgf

% Plastic Centroid
PC = (0.85*b*h^2*fc/2 + sum(as.*d*(fy - fc)))/(P0); % cm

% Neutral Axis depth (c)
beta1_val = beta1(fc);
syms c_
assume(c_, 'positive');
es = (c_ - d)/c_*ecu;
eqn = 0.85*fc*beta1_val*c_*b + sum(as.*sigmaReinf(es,fy,Es)) - Pu;
c = double(solve(eqn,c_)); % cm

% Steel deformation
es = (c - d)/c*ecu; % -

% Forces
sigma_Reinf = sigmaReinf(es,fy,Es); % kgf/cm^2                              % Steel reinforce bars layers stress
f_steel = as.*sigma_Reinf; % kgf                                            % Steel reinforce bars layers forces
Cc = 0.85*fc*beta1_val*c*b; % kgf                                           % Concrete force

% Flexural Strength
Mn = Cc*(PC - beta1_val*c/2) + sum((PC - d).*f_steel); % kgf-cm             % Nominal flexural strength
phi_val = phi(min(es));                                                     % Strength reduction factor
phiMn = phi_val*Mn; % kgf-cm                                                % Design flexural strength

%% Display results
% Now in tonf
Mn_tonf_m = Mn/1000/100;
phiMn_tonf_m = phiMn/1000/100;
Cc_tonf = Cc/1000;
f_steel_tonf = f_steel/1000;
% Display results
for i = 1:length(nBars)
    fprintf('\nRefuerzo %.0f: %.0fphi%.0f a %.0f',i,nBars(i),diams(i)*10,d(i))
end
fprintf('\n c = %.2f [cm]\n \n', c)
tabla = table();
tabla.d = d;
tabla.as = as;
tabla.es = es;
tabla.sigma_steel = sigma_Reinf;
tabla.f_steel = f_steel;
disp(tabla)
fprintf('d [cm] | as [cm2]| es [-] | sigma_steel [kgf/cm2] | f_steel [tonf]\n')
fprintf('Cc = %.2f [tonf]\n',Cc/1000)
fprintf('Cc + f_steels- Pu = %.2f [tonf]\n\n', Cc/1000 + sum(f_steel)/1000 - Pu);
tabla = table();
tabla.Mn_tonf_m = Mn/1000/100;
tabla.Phi = phi_val;
tabla.phiMn_tonf_m = phiMn/1000/100;
disp(tabla)

%% Functions
function sigma = sigmaReinf(es, fy, Es)
    if isa(es,'sym')
        sigma = sym(zeros(length(es),1));
        for i = 1:length(es)
        sigma(i) = piecewise(abs(es(i)) < fy/Es, es(i)*Es, abs(es(i)) >= fy/Es, sign(es(i))*fy);
        end
    elseif isa(es,'double')
        sigma = zeros(length(es),1);
        condition = abs(es) < fy/Es;
        sigma(condition) = es(condition)*Es;
        sigma(~condition) = sign(es(~condition))*fy;
    end
end

function beta1_val = beta1(fc)
    fcMPa = fc/10;
    if fcMPa >= 17 && fcMPa <= 28
        beta1_val = 0.85;
    elseif fcMPa > 28 && fcMPa < 55
        beta1_val = 0.85 - (0.05/7)*(fcMPa-28);
    else
        beta1_val = 0.65;
    end
end

function phi = phi(es)
    if abs(es) < 0.002
        phi = 0.65;
    elseif abs(es) > 0.005
        phi = 0.9;
    else
        phi = 0.65 + (abs(es)-0.002)*(250/3);
    end
end
