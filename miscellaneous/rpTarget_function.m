function residue = rpTarget_function(rp, vInf_minus, vInf_plus, muPlanet, deltaReq)
%
% DESCRIPTION:
%   Evaluates the residual between the required deflection angle (deltaReq) 
%   and the natural gravitational deflection achievable at a specific 
%   pericenter radius (rp).
%   This function is intended to be passed to a zero-finding algorithm 
%   (e.g., fzero) to determine the exact ballistic flyby radius.
%
% INPUT:
%   rp[1]           Radius of pericenter (variable to optimize) [km]
%   vInf_minus[3x1] Incoming hyperbolic excess velocity vector [km/s]
%   vInf_plus[3x1]  Outgoing hyperbolic excess velocity vector [km/s]
%   muPlanet[1]     Gravitational parameter of the flyby planet [km^3/s^2]
%   deltaReq[1]     Required total deflection angle [rad]
%
% OUTPUT:
%   residue[1]      Difference between natural and required deflection [rad]
%                   (residue = 0 implies the solution is found)
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------
    vInf_m = norm(vInf_minus);
    vInf_p = norm(vInf_plus);

    e_m = 1 + (rp*vInf_m^2) / muPlanet;
    e_p = 1 + (rp*vInf_p^2) / muPlanet;

    delta_m = 2*asin(1/e_m);
    delta_p = 2*asin(1/e_p);

    deltaResult = 0.5*delta_m + 0.5*delta_p;

    residue = deltaResult - deltaReq;
end