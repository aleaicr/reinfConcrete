function [l0_dispuesto,n] = l0_dispuestoFunc(l0_norma, s1_disp)
    % Largo mínimo para el primer estribo (5 cm desde el principio)
    l0_estribos = 0;
    
    % Inicialmente, asumimos un solo estribo
    n = 0;
    
    % Mientras el largo calculado no sea mayor al mínimo
    while l0_estribos <= l0_norma
        % Incrementa "n" para verificar si el largo es suficiente
        n = n + 1;
        
        % Calcula el largo considerando "n" estribos
        l0_estribos = 5 + (n - 1)*s1_disp;
    end
    n = n-1;
    % El resultado será el largo calculado
    l0_dispuesto = l0_estribos;
end
