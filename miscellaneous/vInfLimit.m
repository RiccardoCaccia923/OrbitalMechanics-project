function [vInf_dep_max,vInf_ga_max,vInf_arr_max] = vInfLimit(data,constraints)
%
% DESCRIPTION:
%   Establishes the maximum allowable hyperbolic excess velocities (v_inf) 
%   for the departure, gravity assist, and arrival phases based on physical 
%   and operational constraints.
%   - Departure: Derived from the maximum launcher characteristic energy (C3)
%     plus a safety margin.
%   - Flyby: Derived from the maximum physical turning capability
%     at the minimum safe altitude.
%   - Arrival: Set to a fixed operational upper limit for rendezvous/capture.
%
% INPUT:
%   data[struct]        Structure containing planetary physical constants:
%                       .r_Mars     [km]       Mars mean radius
%                       .h_atmMars  [km]       Atmosphere height
%                       .muMars     [km^3/s^2] Mars gravitational parameter
%   constraints[struct] Structure containing mission constraints:
%                       .rocketC3   [km^2/s^2] Max launcher energy
%                       .deltaMin   [rad]      Min required deflection angle
%
% OUTPUT:
%   vInf_dep_max[1]     Maximum departure velocity limit [km/s]
%   vInf_ga_max[1]      Maximum flyby incoming velocity limit [km/s]
%   vInf_Arr_max[1]     Maximum arrival velocity limit [km/s]
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------
    % launcher limit
    safeMarginFactor1 = 1.2;                                         
    vInf_dep_max = sqrt(constraints.rocketC3)*safeMarginFactor1;                   
    
    % turning angle limit
    safeMarginFactor2 = 50;                                     
    rpMin = data.r_Mars + data.h_atmMars + safeMarginFactor2;                       
    vInf_ga_max = sqrt((data.muMars/rpMin)*(1/(sin(constraints.deltaMin/2))-1));   

    % rendezvous limit
    vInf_arr_max = 20;

    fprintf('\n--- vInf Limitations for Leg 1---\n');
    fprintf('vInf max @Jupiter: %.3f km/s \n',vInf_dep_max);
    fprintf('vInf max @Mars: %.3f km/s \n',vInf_ga_max);
    fprintf('vInf max @Asteroid: %.3f km/s \n',vInf_arr_max);
end
