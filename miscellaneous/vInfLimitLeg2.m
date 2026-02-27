function [vInf_ga_max,vInf_arr_max] = vInfLimitLeg2(data,constraints)
    
    % turning angle limit
    safeMarginFactor2 = 50;                                     
    rpMin = data.r_Mars + data.h_atmMars + safeMarginFactor2;                      %[km] 
    vInf_ga_max = sqrt((data.muMars/rpMin)*(1/(sin(constraints.deltaMin/2))-1));   %[km/s]

    % rendezvous limit
    vInf_arr_max = 20;

    fprintf('\n--- vInf Limitations for Leg2 ---\n');
    fprintf('vInf max @gravityAssistPlanet: %.3f km/s \n',vInf_ga_max);
    fprintf('vInf max @arrivalAsteroid: %.3f km/s \n',vInf_arr_max);
end

% !!! -> GO TO vInfLimit