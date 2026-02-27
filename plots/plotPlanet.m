function plotPlanet(planetTexture,planetRadius,planetCenter)
%
% DESCRIPTION:
%   Generates a high-resolution 3D spherical surface representing a planet.
%
% INPUT:
%   planetTexture   [char]   Filename or path of the texture image (e.g., 'mars.jpg')
%   planetRadius    [1]      Radius of the sphere (in plot units, e.g., km or AU)
%   planetCenter    [3x1]    (Optional) Center coordinates [x, y, z].
%                            Default is [0, 0, 0] if omitted or empty.
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------
if nargin < 3 || isempty(planetCenter)
        planetCenter = [0, 0, 0];
end

% Generate the sphere
[xP,yP,zP] = sphere(200);   
   
% Load image
Img = imread(planetTexture);
Img = flipud(Img);      % fix nord axis as z positive

X = xP * planetRadius + planetCenter(1);
Y = yP * planetRadius + planetCenter(2);
Z = zP * planetRadius + planetCenter(3);

% Create surface with Earth texture
surface(X, Y, Z, ...
            'CData', Img, ...
            'FaceColor', 'texturemap', ...
            'EdgeColor', 'none', ...
            'FaceLighting', 'none');

end