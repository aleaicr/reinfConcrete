function phi = phi(es)
    if abs(es) < 0.002
        phi = 0.65;
    elseif abs(es) > 0.005
        phi = 0.9;
    else
        phi = 0.65 + (abs(es)-0.002)*(250/3);
    end
end