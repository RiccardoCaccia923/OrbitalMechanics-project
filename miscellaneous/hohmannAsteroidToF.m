function [TOF, TOF_days] = hohmannAsteroidToF(a_Dep,a_Arr)
% hohmannAsteroidToF Computes the duration of a transfer between two semi-major axes.
%
%
% DESCRIPTION:
%   Calculates strictly the Time of Flight (half-period) of a transfer orbit 
%   connecting two distances (or semi-major axes). 
%   Unlike the full Hohmann function, this does not compute Delta V, making 
%   it useful for quick timing estimates or phasing analysis where orbital 
%   velocities are not immediately required.
%   IT HAS BEEN USED ONLY IN A PRELIMINARY PHASE TO UNERSTAND APPROXIMATELY
%   THE TOF
%
% INPUT:
%   a_Dep[1]    Semi-major axis (or radius) of departure [km]
%   a_Arr[1]    Semi-major axis (or radius) of arrival [km]
%
% OUTPUT:
%   TOF[1]      Time of Flight [s]
%   TOF_days[1] Time of Flight [days]
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------
    mu = astroConstants(4);

    a_t = (a_Dep + a_Arr) / 2;

    TOF = pi * sqrt(a_t^3 / mu); % Tempo in secondi
    TOF_days = TOF / 86400;

    fprintf('\n--- Hohmann transfer from Mars to Asteroid ---\n');
    fprintf('Hohmann Time of Flight : %.1f days (%.1f years)\n', TOF_days, TOF_days/365.25);

end