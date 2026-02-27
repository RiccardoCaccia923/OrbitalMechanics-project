function plotPlanetUnit(plColor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotPlanetUnit: this function creates a sphere covered with generic planet layout 
%           in order to complete orbit plots
% 
%                      GO TO plotPlanet 
% CONTRIBUTORS:
%   Riccardo Caccia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the sphere
[xP,yP,zP] = sphere(50);       

planet_surf = surf(xP, yP, zP); 

set(planet_surf, 'FaceColor', plColor, ... 
                 'EdgeColor', 'none', ...       
                 'FaceAlpha', 0.6, ...          
                 'DiffuseStrength', 1, ...
                 'SpecularStrength', 0.2, ...
                 'HandleVisibility', 'off');

light('Position', [10 10 10], 'Style', 'infinite');
lighting gouraud;

end