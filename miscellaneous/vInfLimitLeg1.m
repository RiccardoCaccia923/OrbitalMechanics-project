function [vInf_dep_max,vInf_ga_max] = vInfLimitLeg1(data,constraints)
    
    % launcher limit
    safeMarginFactor1 = 1.2;                                         
    vInf_dep_max = sqrt(constraints.rocketC3)*safeMarginFactor1;                   %[km/s]
    
    % turning angle limit
    safeMarginFactor2 = 50;                                     
    rpMin = data.r_Mars + data.h_atmMars + safeMarginFactor2;                      %[km] 
    vInf_ga_max = sqrt((data.muMars/rpMin)*(1/(sin(constraints.deltaMin/2))-1));   %[km/s]

    fprintf('\n--- vInf Limitations for Leg 1---\n');
    fprintf('vInf max @departurePlanet: %.3f km/s \n',vInf_dep_max);
    fprintf('vInf max @gravityAssistPlanet: %.3f km/s \n',vInf_ga_max);
end

% !!! -> GO TO vInfLimit