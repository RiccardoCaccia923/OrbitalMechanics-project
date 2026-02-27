function kep = car2kep(mu, s)
%
% DESCRIPTION
%   Converts Cartesian state vector [r; v] -> Keplerian elements
%
% INPUT
%   mu  : gravitational parameter [km^3/s^2]
%   s   : column vector [6x1] = [rx ry rz vx vy vz]^T
%
% OUTPUT (struct)
%   kep.a     : semi-major axis (negative if hyperbolic)
%   kep.e     : eccentricity (magnitude)
%   kep.i     : inclination [rad]
%   kep.raan  : RAAN Ω [rad]
%   kep.w     : argument of pericenter ω [rad]
%   kep.theta : true anomaly θ [rad]
%
% Notes on singularities/edge cases:
% - If e ~ 0 (circular), ω is undefined -> set to 0, and θ becomes
%   the argument of latitude u (or true longitude for equatorial orbits).
% - If i ~ 0 (equatorial), RAAN is undefined -> set to 0, and ω is
%   measured from the x-axis to the pericenter axis in the orbital plane.
%
% CONTRIBUTORS:
%   Riccardo Caccia
%--------------------------------------------------------------------------


s = s(:);
r = s(1:3); v = s(4:6);

rNorm = norm(r);
vNorm = norm(v);
if rNorm < eps
    error('Norma di r nulla: stato non valido.');
end

khat = [0;0;1];
tol_i = 1e-12;         % tolleranza per inclinazione/RAAN
tol_e = 1e-12;         % tolleranza per eccentricità (circolare)
tol_vr = 1e-12;        % tolleranza per velocità radiale (perigeo/apogeo)

% -------------------- grandezze ausiliarie --------------------
h = cross(r, v);
hNorm = norm(h);
if hNorm < 1e-14
    error('Momento angolare ~ 0: traiettoria rettilinea non orbitale.');
end
uh = h / hNorm;

ur = r / rNorm;                 % versore radiale
vt = v - dot(v,ur)*ur;          % componente trasversale vettoriale
ut = vt / max(norm(vt), realmin); % versore trasversale (se v ~ radiale evita NaN)
vr = dot(v, ur);                % velocità radiale

% vettore nodale
N = cross(khat, h);
Nnorm = norm(N);
if Nnorm > tol_i
    uN = N / Nnorm;
else
    uN = [1;0;0]; % segnaposto per casi quasi equatoriali
end

% eccentricità (vettore e)
% e_vec = ( (v^2 - mu/r)*r - (r·v)*v ) / mu       (formula più stabile)
e_vec = ((vNorm^2 - mu/rNorm)*r - dot(r,v)*v) / mu;
eNorm = norm(e_vec);

% -------------------- elementi classici --------------------
% inclinazione i
ci = dot(uh, khat);        % cos(i)
ci = max(-1, min(1, ci));
i = acos(ci);

% RAAN Ω = atan2(N_y, N_x) se non equatoriale
if i > tol_i
    raan = atan2(uN(2), uN(1));
    raan = mod(raan, 2*pi);
else
    raan = 0; % indeterminato in equatoriale
end

% argomento del pericentro ω e anomalia vera θ
if eNorm > tol_e
    ue = e_vec / eNorm;
    % ω misurato lungo la linea dei nodi: usa atan2 per il quadrante
    if i > tol_i
        % ω = angle(uN -> ue) nel piano orbitale
        sinw = dot(cross(uN, ue), uh);
        cosw = dot(uN, ue);
        w = atan2(sinw, cosw);
    else
        % equatoriale: line of nodes non definita; usa asse x come riferimento
        ex = [1;0;0];
        sinw = dot(cross(ex, ue), uh);
        cosw = dot(ex, ue);
        w = atan2(sinw, cosw);
    end
    w = mod(w, 2*pi);

    % θ = angle(ue -> ur) nel piano orbitale (atan2 per quadrante)
    sinth = dot(cross(ue, ur), uh);
    costh = dot(ue, ur);
    theta = atan2(sinth, costh);
    theta = mod(theta, 2*pi);

    % opzionale: forzatura per perigeo/apogeo se vr ~ 0
    if abs(vr) < tol_vr
        if costh >= 0
            theta = 0;        % perigeo
        else
            theta = pi;       % apogeo
        end
    end
else
    % orbita (quasi) circolare -> ω indeterminato
    w = 0;

    if i > tol_i
        % argomento di latitudine u = angle(uN -> ur)
        sinu = dot(cross(uN, ur), uh);
        cosu = dot(uN, ur);
        theta = atan2(sinu, cosu);  % qui theta rappresenta u
    else
        % equatoriale + circolare: true longitude dall’asse x
        theta = atan2(ur(2), ur(1));
    end
    theta = mod(theta, 2*pi);
end

% semiasse maggiore a (valido per tutte le coniche: ellittica, iperbolica)
% vis-viva: a = 1 / ( 2/r - v^2/mu ); per parabola (e=1) -> a = +Inf
a = 1 / ( 2/rNorm - (vNorm^2)/mu );

% -------------------- output --------------------
kep.a     = a;
kep.e     = norm(e_vec);
kep.i     = i;
kep.raan  = raan;
kep.w     = w;
kep.theta = theta;
end
