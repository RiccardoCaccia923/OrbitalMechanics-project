function time = hyperbolicTimeLaw(v_inf, r_p, mu, r_lim)
%
% DESCRIPTION:
%   Calculates the time required for a spacecraft on a hyperbolic orbit 
%   to travel from the pericenter (r_p) to a specific radial distance
%   (r_lim).
%
% INPUT:
%   v_inf[1]        Hyperbolic excess velocity [km/s]
%   r_p[1]          Radius of pericenter [km]
%   mu[1]           Gravitational parameter of the central body [km^3/s^2]
%   r_lim[1]        Limit radius (e.g., Radius of SOI) [km]
%
% OUTPUT:
%   time[1]         Time of flight from pericenter to r_lim [s]
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------
    a = -mu / v_inf^2;        
    e = 1 + (r_p * v_inf^2)/mu; 
    
    cosh_F = (1 - r_lim/a) / e;
    if cosh_F < 1
        cosh_F = 1; 
    end
    
    F = acosh(cosh_F);

    M = e * sinh(F) - F;
    
    n = sqrt(mu / abs(a)^3);
    
    time = M / n;
end