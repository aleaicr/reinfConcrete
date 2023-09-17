function resultado = redondearEsp(numero)
 par_cercano = 2 * floor(numero / 2);
    multiplo5_cercano = 5 * floor(numero / 5);
    
    % Calcula las diferencias absolutas
    diferencia_par = abs(numero - par_cercano);
    diferencia_multiplo5 = abs(numero - multiplo5_cercano);
    
    % Compara las diferencias y elige el redondeo más cercano hacia abajo
    if diferencia_par < diferencia_multiplo5
        resultado = par_cercano;
    else
        resultado = multiplo5_cercano;
    end
end
