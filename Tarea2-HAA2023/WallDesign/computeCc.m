function [Cc, Cc_centroid] = computeCc(b, h, fc, c, beta1_val)
% This function computes the contribution of the concrete to the flexural
% strength (Cc) and the centroid of the distributed force (Cc_centroid)
%
% INPTUS
% b: vector of widths of the section
% h: vector of heights of the section
% c: neutral axis distance (to the top of the section)
% beta1_val: value of beta1
%
% OUTPUTS
% Cc:
% Cc_centroid:
%
% NOTES
%
%

% Previous
h_ = cumsum(h); 

%% If c is double
% Determinar rango de c
h_ = cumsum(h); 
h_inf = max(h_(h_<=beta1_val*c));
h_sup = min(h_(h_>=beta1_val*c));

% Determinar todas las areas que contribuyen a Cc
area = [b(h_<=h_inf).*h(h_<= h_inf); (beta1_val*c-h_(h_==h_inf))*b(h_==h_sup)];

% Determinar Cc y centroide de Cc
Cc = 0.85*fc*area; % kgf (vector)
aux = h;
aux(end) = [];
aux = [0; cumsum(aux)];
h_centroids = h/2 + aux;
Cc_centroid = sum(area.*h_centroids)/sum(area);

%% If c is symbolic
h_length = length(h_);
Cc = syms(zeros(h_length,1));


end

