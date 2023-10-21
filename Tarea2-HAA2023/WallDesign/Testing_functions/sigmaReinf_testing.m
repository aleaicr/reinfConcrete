function sigma = sigmaReinf_testing(es, fy, Es)
% Testing if I can get sigma for a matriz es
% [m,n] = size(es)
% m = number of rows, number of reinforcement layers
% n = number of different 'es' I want to try (eg.: Moment-Curvature Diagram)
% this shouldn't work because piecewise does not accept matrix

if isa(es,'sym')
    sigma = piecewise(abs(es) < fy/Es, es*Es, abs(es) >= fy/Es, sign(es)*fy); % Symbolic vector of symbolic piecewises of a symbolic variable
elseif isa(es,'double')
    sigma = zeros(size(es));
    condition = abs(es) < fy/Es;
    sigma(condition) = es(condition)*Es;
    sigma(~condition) = sign(es(~condition))*fy;
end
end
