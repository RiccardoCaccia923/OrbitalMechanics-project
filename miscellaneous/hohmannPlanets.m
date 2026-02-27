function [Dv, TOF, TOF_days] = hohmannPlanets(R_Dep,R_Arr)
% hohmannPlanets Calculates cost and duration of a Hohmann transfer between circular orbits.
%
% DESCRIPTION:
%   Computes the total Delta V and Time of Flight for a bi-elliptic Hohmann 
%   transfer between two concentric, coplanar, circular orbits.
%   It calculates the impulsive maneuvers required at departure and arrival 
%   to shift from the initial radius (R_Dep) to the final radius (R_Arr).
%   IT HAS BEEN USED ONLY IN A PRELIMINARY PHASE TO UNERSTAND APPROXIMATELY
%   THE TOF
%
% INPUT:
%   R_Dep[1]    Radius of the departure circular orbit [km]
%   R_Arr[1]    Radius of the arrival circular orbit [km]
%
% OUTPUT:
%   Dv[1]       Total Delta V cost (Departure burn + Arrival burn) [km/s]
%   TOF[1]      Time of Flight (transfer duration) [s]
%   TOF_days[1] Time of Flight [days]
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------
    mu = astroConstants(4); % km^3/s^2
    
    V_JupCirc  = sqrt(mu / R_Dep);
    V_MarsCirc = sqrt(mu / R_Arr);
    
    a_t = (R_Dep + R_Arr) / 2;
    
    V_tDep = sqrt(mu * (2/R_Dep - 1/a_t));  
    V_tArr = sqrt(mu * (2/R_Arr - 1/a_t)); 
    
    Dv1 = abs(V_tDep - V_JupCirc);  
    Dv2 = abs(V_tArr - V_MarsCirc);
    
    Dv = Dv1 + Dv2;
    
    TOF = pi * sqrt(a_t^3 / mu);
    TOF_days = TOF / 86400;
    
    fprintf('\n--- Hohmann transfer from Jupiter to Mars ---\n');
    fprintf('Hohmann DV : %.3f km/s \n', Dv);
    fprintf('Hohmann Time of Flight : %.1f days (%.1f years)\n', TOF_days, TOF_days/365.25);
end