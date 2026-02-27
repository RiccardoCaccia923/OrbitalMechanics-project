function [Tsyn,Tsyn_days,Tsyn_years] = computeTsynLeg1(a1, a2, mu)
% computeTsynLeg1 Computes the synodic period between two orbiting bodies.
%
% DESCRIPTION:
%   Calculates the synodic period of bodies in the first interplanetary leg
% INPUT:
%   a1[1]   Semi-major axis of the first body [km]
%   a2[1]   Semi-major axis of the second body [km]
%   mu[1]   Gravitational parameter of the central body [km^3/s^2]
%
% OUTPUT:
%   Tsyn[1]         Synodic period [s]
%   Tsyn_days[1]    Synodic period [days] (Tsyn / 86400)
%   Tsyn_years[1]   Synodic period [years] (Tsyn_days / 365.25)
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------

    T1 = 2 * pi * sqrt(a1^3 / mu);
    T2 = 2 * pi * sqrt(a2^3 / mu);

    Tsyn = (T1 * T2) / abs(T1 - T2);

    Tsyn_days = Tsyn / 86400;       
    Tsyn_years = Tsyn_days / 365.25; 

    fprintf('\n--- Synodic Period ---\n');
    fprintf('Jupiter-Mars Synodic Period: %.1f days (%.2f years)\n', ...
        Tsyn_days, Tsyn_years); % change if planets change
end