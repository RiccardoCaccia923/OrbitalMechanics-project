function s = kep2car(mu, kep)
% DESCRIPTION
%   Converts Keplerian elements -> Cartesian state [r; v]
%
% INPUT:
%   mu        : gravitational parameter [km^3/s^2]
%   kep       : struct with fields
%               .a     semi-major axis [km]
%               .e     eccentricity [-]
%               .i     inclination [rad]
%               .raan  right ascension of the ascending node Ω [rad]
%               .w     argument of pericenter ω [rad]
%               .theta true anomaly θ [rad]
%
% OUTPUT:
%   s         : [6x1] Cartesian state = [r; v]
%               r = [rx; ry; rz] [km]
%               v = [vx; vy; vz] [km/s]
%
% Valid for elliptical and hyperbolic orbits.
% Automatically handles equatorial / circular cases (robust/singularity-free).
%
% CONTRIBUTORS:
%   Riccardo Caccia
% -------------------------------------------------------------------------

a     = kep.a;
e     = kep.e;
i     = kep.i;
raan  = kep.raan;
w     = kep.w;
theta = kep.theta;

% ---------------- controlli input -----------------
if e < 0
    error('Eccentricità negativa non valida.');
end
if abs(i) > pi
    error('Inclinazione fuori intervallo [-pi, pi].');
end

% ---------------- coordinate nel piano orbitale ----------------
% r e v nel sistema PQW (perifocale)
p = a*(1 - norm(e)^2);           % parametro orbitale [km]
r_pf = (p / (1 + norm(e)*cos(theta))) * [cos(theta); sin(theta); 0];

v_pf = sqrt(mu/p) * [-sin(theta);
                      norm(e) + cos(theta);
                      0];

% ---------------- matrici di rotazione ----------------
R3_W = [ cos(raan)  sin(raan)  0;
        -sin(raan)  cos(raan)  0;
              0           0     1 ];

R1_i = [1     0           0;
        0  cos(i)   sin(i);
        0 -sin(i)   cos(i)];

R3_w = [ cos(w)  sin(w)  0;
        -sin(w)  cos(w)  0;
              0       0  1 ];

% Rotazione complessiva PQW -> IJK (inerziale)
Q_pX = (R3_W') * (R1_i') * (R3_w');
% oppure più esplicitamente:
% Q_pX = R3(-raan) * R1(-i) * R3(-w)

% ---------------- vettori in sistema inerziale ----------------
r = Q_pX * r_pf;
v = Q_pX * v_pf;

r = snapzero(r);
v = snapzero(v);

% ---------------- output ----------------
s = [r; v];
end
