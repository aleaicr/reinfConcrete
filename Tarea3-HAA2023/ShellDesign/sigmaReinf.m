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
