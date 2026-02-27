function plotGAOpt(data,GA)
% plotGA Visualizes the planetocentric hyperbolic trajectory of the gravity assist.
%
%
% DESCRIPTION:
%   Generates a detailed 3D plot of the flyby maneuver in the Mars-centered
%   frame (Planetocentric).
%   The function:
%   1. Reconstructs the incoming and outgoing hyperbolic arcs starting from 
%      the pericenter conditions derived in the optimization phase.
%   2. Applies the necessary rotation to the incoming velocity vector to 
%      account for the powered maneuver (turn angle correction).
%   3. Visualizes the planet (with texture), the trajectory arcs, the 
%      pericenter location, and the theoretical asymptotes.
%   4. Prints the orbital elements of both hyperbolas to the 
%      console for verification.
%
% INPUT:
%   GA[struct]      Gravity Assist solution structure (output of computeGA):
%                   .rp, .vp_minus, .vp_plus, .vInf_minus, .vInf_plus, 
%                   .delta_maneuver (optional)
%   data[struct]    Planetary data structure:
%                   .muMars [km^3/s^2], .r_Mars [km], texture filename
%
% OUTPUT:
%   The function prints the calculated hyperbolic parameters:
%   - Incoming Hyperbola: Semi-major axis (a-), Eccentricity (e-)
%   - Outgoing Hyperbola: Semi-major axis (a+), Eccentricity (e+)
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------
rp = GA.rp;
vp_minus = GA.vp_minus;
vp_plus = GA.vp_plus;
vInf_minus = GA.vInf_minus;
vInf_plus = GA.vInf_plus;

r0_perif = [rp; 0; 0];
v0_IN  = [0; vp_minus; 0];  
v0_OUT = [0; vp_plus; 0];

if isfield(GA, 'delta_maneuver')
    turn_angle = GA.delta_maneuver; 
else
    turn_angle = 0;
end
u_rot = [0; 0; 1]; 
v0_IN = rotateVector(v0_IN, u_rot, -turn_angle);

sp_OUT = [r0_perif;v0_OUT];
tspanfw = linspace(0,1e4,1e5);
[~,sfw] = twoBodyProblemNoPlotHyp(data.muMars,sp_OUT,[],tspanfw);
sp_IN = [r0_perif;v0_IN];
tspanbw = linspace(0,-1e4,1e5);
[~,sbw] = twoBodyProblemNoPlotHyp(data.muMars,sp_IN,[],tspanbw);

% plot hyperbolic arcs
scale = data.r_Mars;
% label_unit = 'R_p';
% plColor = [0.7 0.7 0.7]; 
% uvpl = Vplanet/norm(Vplanet);

figure('FileName','Planetocentric Plot')
plot3(sbw(:,1)/scale,sbw(:,2)/scale,sbw(:,3)/scale,'LineWidth',2,'color','[0.47, 0.67, 0.19]')
hold on
plot3(sfw(:,1)/scale,sfw(:,2)/scale,sfw(:,3)/scale,'LineWidth',2,'color','[0.4940, 0.1840, 0.5560]')
plot3(r0_perif(1)/scale,r0_perif(2)/scale,r0_perif(3)/scale,'bo','LineWidth',1.5)
plotPlanet('MarsTexture.jpg',data.r_Mars/scale)
plot3([-5,5],[0,0],[0,0],'k--','HandleVisibility', 'off')
% quiver3(0,0,0,uvpl(1),uvpl(2),uvpl(3),2,'k','LineWidth',1.5)
grid on
axis equal
title('Flyby Hyperbola (Planetocentric Frame)')
xlim([-2,2.5])
ylim([-4,4])
zlim([-3,3])
xlabel('x [Rpl]'); 
ylabel('y [Rpl]'); 
zlabel('z [Rpl]');
legend('Incoming Arc','Outgoing Arc','Pericenter','Mars','Location', 'best');

% compute asintotes
vInf_m = norm(vInf_minus);
plot_limit_x = -10*scale ;
a_min = -data.muMars / vInf_m^2;
e_min = 1 + (rp * vInf_m^2) / data.muMars;
Cx_min = a_min * e_min;
nu_inf_min = acos(-1/e_min);
m_in = tan(nu_inf_min);

x_raw = linspace(plot_limit_x, -Cx_min, 100);
y_raw = -m_in * (x_raw + Cx_min); 
z_raw = zeros(size(x_raw)); 

R = [cos(-turn_angle), -sin(-turn_angle); 
     sin(-turn_angle),  cos(-turn_angle)];
pts_raw = [x_raw; y_raw]; 
pts_rot = R * pts_raw;    

x_asym_in = pts_rot(1,:);
y_asym_in = pts_rot(2,:);

plot(x_asym_in/scale, y_asym_in/scale, '--','color','[0.47, 0.67, 0.19]', ...
     'LineWidth', 1.5, 'DisplayName', 'Incoming Hyp. Asymptote');

vInf_p = norm(vInf_plus);
plot_limit_x = -10*scale;
a_plus = -data.muMars / vInf_p^2;
e_plus = 1 + (rp * vInf_p^2) / data.muMars;
Cx_plus = a_plus * e_plus; 
nu_inf_plus = acos(-1/e_plus);
m_out = tan(-nu_inf_plus);
x_asym_out = linspace(plot_limit_x, -Cx_plus, 100);
y_asym_out = -m_out * (x_asym_out + Cx_plus);
plot(x_asym_out/scale, y_asym_out/scale, '--','color','[0.4940, 0.1840, 0.5560]', ...
     'LineWidth', 1.5, 'DisplayName', 'Outgoing Hyp. Asymptote');

legend('Location', 'best');

fprintf('\n---ORBITAL PARAMETERS---\n');

fprintf('INcoming hyperbole:\n');
fprintf('  > Semi-Major Axis (a-): %.4f km\n', a_min);
fprintf('  > Eccentricity (e-):      %.4f\n', e_min);
fprintf('\n');

fprintf('OUTgoing hyperbole:\n');
fprintf('  > Semi-MAjor Axis (a+): %.4f km\n', a_plus);
fprintf('  > Eccentricity (e+):      %.4f\n', e_plus);
fprintf('------------------------------------------------\n');

end