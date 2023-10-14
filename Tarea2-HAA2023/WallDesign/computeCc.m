function [Cc, Cc_centroid] = computeCc(b, h, fc, c, beta1_val)
% Tarea 2 - Hormigón Armado Avanzado
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras R. - Gabriel Ramos V.
%%
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
% This function solves for symbolic 'c' or for a given value of 'c'
%

%% Previous Calcs
h_ = cumsum(h); 
beta1c = beta1_val*c;                                                           % puede ser double o simbólico

%% If c is symbolic (para getMn --> hay que calcular 'c')
if isa(c,'sym')
    h_length = length(h_);
    area = sym(zeros(h_length,1));
    h_inf = 0;
    for i = 1:h_length
        h_sup = h_(i);
        area(i) = piecewise(beta1c <= h_inf, 0, beta1c > h_inf & beta1c < h_sup, (beta1c-h_inf)*b(i), beta1c > h_sup, (h_sup-h_inf)*b(i));
        h_inf = h_sup;
    end
    Cc = 0.85*fc*area; % kgf (vector)
%     aux = h;
%     aux(end) = [];
%     aux = [0; cumsum(aux)];
%     h_centroids = h/2 + aux;
%     Cc_centroid = sum(area.*h_centroids)/sum(area);
    Cc_centroid = 0;
end

%% If c is double  (para getMn_esBased --> ya conozco 'c')
if isa(c,'double')
    % Determinar rango de c
    h_inf = max(h_(h_ <= beta1c));
    h_sup = min(h_(h_ >= beta1c));
    % Determinar todas las areas que contribuyen a Cc
    if isempty(h_inf)       % Si estoy en la primera zona
        area = beta1c*b(1);
    else
        area = [h(h_<=h_inf).*b(h_<= h_inf); (beta1c-h_(h_==h_inf)).*b(h_==h_sup)];
    end
    % Determinar Cc y centroide de Cc
    Cc = 0.85*fc*area; % kgf (vector)
    aux = h;
    aux(end) = [];
    aux = [0; cumsum(aux)];
    h_centroids = h/2 + aux;
    Cc_centroid = sum(area.*h_centroids)/sum(area);
end

end

