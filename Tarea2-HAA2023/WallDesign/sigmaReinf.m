function sigma = sigmaReinf(es, fy, Es)
% Computes the stress of each layer of steel
%
% Inputs
% es: (vector) strain of each layer of steel
% fy: Steel's strength
% Es: Steel's modulus of elasticity
%
% Outputs:
% sigma: Stress of each layer of steel
%
% Notes
% * 'es' can be symbolic vector of an expression in function of 'c' 
% (c is the depth of the neutral axis). This, for use in the force
% equilibrium of the section to determine the value of 'c'.
% * 'es' can be a double vector 
if isa(es,'sym')
    sigma = sym(zeros(length(es),1));                                   % Symbolic vector of symbolic piecewises of a symbolic variable
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
