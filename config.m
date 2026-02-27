function [constants,data,constraints,Leg1,Leg2] = config()
% config Initializes the simulation environment and mission parameters.
%
% DESCRIPTION:
%   Acts as the centralized configuration loader for the mission analysis.
%   This function:
%   1. Defines universal physical constants (using astroConstants).
%   2. Sets planetary and target data (IDs, physical properties).
%   3. Establishes operational constraints (Launcher C3, Min Flyby Angle).
%
% INPUT:
%   None.
%
% OUTPUT:
%   constants[struct]   Physical constants (muSun, AU).
%   data[struct]        Body-specific data (IDs, Radii, mu, Atmospheres).
%   constraints[struct] Operational constraints (rocketC3, deltaMin).
%   Leg1[struct]        Velocity upper bounds for the first leg.
%   Leg2[struct]        Velocity upper bounds for the second leg.
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
%--------------------------------------------------------------------------

% CONSTANTS 
constants.muSun = astroConstants(4);
constants.AU = astroConstants(2);

% DATA 
data.r_Sun = 696340;
data.muSun = astroConstants(4);
data.idDep = 5;
data.idGA = 4;
data.asteroidNumber = 448628;
data.R_Jup  = 5.204 * constants.AU;  
data.r_Jup = 71492;
data.R_Mars = 1.524 * constants.AU;  
data.muMars = astroConstants(14);
data.r_Mars = 3396;
data.h_atmMars = 120;

% CONSTRAINTS 
constraints.rocketC3 = 100;                                                
constraints.deltaMin = deg2rad(1); 
[Leg1.vInf_dep_max,Leg1.vInf_ga_max,Leg2.vInf_arr_max] = vInfLimit(data,constraints);
Leg2.vInf_ga_max = Leg1.vInf_ga_max;

end