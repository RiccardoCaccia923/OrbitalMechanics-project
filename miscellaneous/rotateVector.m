function vecRot = rotateVector(v,u,delta)
%   DESCRIPTION: 
%    function which uses Rodrigues's formula to rotate a
%    vector v of an angle delta around a unit vector u (clockwise)
%
%   INPUT:
%           v: original vector
%           u: vector to rotate around
%           delta: angle of rotation around u [rad]
%
%   OUTPUT:
%           vecRot: vector result of rotation
%
%   CONTRIBUTORS:
%     Riccardo Caccia
%
% -------------------------------------------------------------------------

%delta = deg2rad(delta);
vecRot = v*cos(delta) + cross(u,v)*sin(delta) + u*dot(u,v)*(1-cos(delta));

end