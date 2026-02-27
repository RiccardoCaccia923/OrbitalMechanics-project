function plot3Opt(constants,data,tDep_opt,tGA_opt,tArr_opt)
% plot3Opt Visualizes the final optimized 3D mission trajectory.
%
% PROTOTYPE:
%   plot3Opt(constants, data, tDep_opt, tGA_opt, tArr_opt)
%
% DESCRIPTION:
%   Generates a comprehensive 3D visualization of the complete mission profile 
%   corresponding to the global optimum found by the solver.
%   The function:
%   1. Computes and plots the full orbits of the Departure body, Flyby body, 
%      and Target Asteroid at the specific optimal epochs.
%   2. Places textured 3D models of the Sun, planets, and a placeholder for 
%      the asteroid at their precise positions.
%   3. Reconstructs the transfer legs (Leg 1 & Leg 2) by solving Lambert's 
%      problem between the optimized timestamps to visualize the actual 
%      spacecraft path.
%
% INPUT:
%   constants[struct]   Physical constants (must include .muSun, .AU)
%   data[struct]        Mission data (IDs, radii, texture filenames)
%   tDep_opt[datetime]  Optimal Departure Date
%   tGA_opt [datetime]  Optimal Gravity Assist Date
%   tArr_opt[datetime]  Optimal Arrival Date
%
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------

mjd2000Dep = date2mjd2000(datetime2vec(tDep_opt));
[kepDep] = uplanet(mjd2000Dep, data.idDep);
kepDep = kepVec2Struct(kepDep);
sDep0 = kep2car(constants.muSun,kepDep);
[~,sDep] = twoBodyProblemNoPlotHyp(constants.muSun,sDep0,1,[]);

mjd2000GA = date2mjd2000(datetime2vec(tGA_opt));
[kepGA] = uplanet(mjd2000GA, data.idGA);
kepGA = kepVec2Struct(kepGA);
sGA0 = kep2car(constants.muSun,kepGA);
[~,sGA] = twoBodyProblemNoPlotHyp(constants.muSun,sGA0,1,[]);

mjd2000Ast = date2mjd2000(datetime2vec(tArr_opt));
kepAst = ephAsteroids( mjd2000Ast, data.asteroidNumber);
kepAst = kepVec2Struct(kepAst);
sAst0 = kep2car(constants.muSun,kepAst);
[~,sAst] = twoBodyProblemNoPlotHyp(constants.muSun,sAst0,1,[]);

figure

scale = constants.AU;

plot3(sDep(:,1)/scale,sDep(:,2)/scale,sDep(:,3)/scale,'LineWidth',1.5)
hold on
plot3(sGA(:,1)/scale,sGA(:,2)/scale,sGA(:,3)/scale,'LineWidth',1.5)
plot3(sAst(:,1)/scale,sAst(:,2)/scale,sAst(:,3)/scale,'color','[0.6, 0.6, 0.6]','LineWidth',1.5)

JupCenter = [sDep(1,1)/scale,sDep(1,2)/scale,sDep(1,3)/scale];
plotPlanet('JupiterTexture.jpg',data.r_Jup*900/scale,JupCenter)

MarsCenter = [sGA(1,1)/scale,sGA(1,2)/scale,sGA(1,3)/scale];
plotPlanet('MarsTexture.jpg',data.r_Mars*8000/scale,MarsCenter)

AstCenter = [sAst(1,1)/scale,sAst(1,2)/scale,sAst(1,3)/scale];
[x_u, y_u, z_u] = sphere(30);
X_S = x_u * 0.1 + AstCenter(1);
Y_S = y_u * 0.1 + AstCenter(2);
Z_S = z_u * 0.1 + AstCenter(3);
surf(X_S, Y_S, Z_S, ...
     'FaceColor', [0.25, 0.25, 0.25], ...
     'EdgeColor', 'none', ...
     'FaceLighting', 'gouraud', ...
     'AmbientStrength', 0.5);

plotPlanet('SunTexture.jpg',data.r_Sun*50/scale)

TofLeg1 = seconds(tGA_opt - tDep_opt);
[~,~,~,~,vLeg1,~,~,~] = lambertMR(sDep0(1:3), sGA0(1:3), TofLeg1, constants.muSun, 0, 0, 0);
sLeg1_0 = [sDep0(1:3);vLeg1'];
[~,sLeg1,~] = twoBodyProblemNoPlotHyp(constants.muSun,sLeg1_0,[],linspace(0,TofLeg1,1e5));
plot3(sLeg1(:,1)/scale,sLeg1(:,2)/scale,sLeg1(:,3)/scale,'color','[0.47, 0.67, 0.19]','LineWidth',1.5)

TofLeg2 = seconds(tArr_opt - tGA_opt);
[~,~,~,~,vLeg2,~,~,~] = lambertMR(sGA0(1:3), sAst0(1:3), TofLeg2, constants.muSun, 0, 0, 0);
sLeg2_0 = [sGA0(1:3);vLeg2'];
[~,sLeg2,~] = twoBodyProblemNoPlotHyp(constants.muSun,sLeg2_0,[],linspace(0,TofLeg2,1e5));
plot3(sLeg2(:,1)/scale,sLeg2(:,2)/scale,sLeg2(:,3)/scale,'color','[0.4940, 0.1840, 0.5560]','LineWidth',1.5)


xlabel('X [AU]')
ylabel('Y [AU]')
zlabel('Z [AU]')
grid on
axis equal
legend('Dep Orbit','GA Orbit','Ast Orbit','Jupiter','Mars','Asteroid','Sun','Leg 1','Leg 2','Location','best')

end