function [kepDep,kepGA,kepAst] = plotAllorbits(constants,data)
% plotAllorbits Visualizes the general orbital configuration of the mission.
%
%
% DESCRIPTION:
%   Generates a 3D plot of the complete orbits for the Departure body, 
%   Gravity Assist body, and Target Asteroid. 
%   The function uses nominal/indicative dates (e.g., 2030 start, 2060 end) 
%   solely to compute and display the geometry, spatial scales, and relative 
%   inclinations of the celestial bodies involved in the mission.
%   This serves as an initial visual check of the system configuration 
%   before starting the trajectory optimization routine.
%
%   FOR BETTER PLOT -> GO TO -> plot3Opt
%
% INPUT:
%   constants[struct]   Physical constants (must include .muSun, .AU)
%   data[struct]        Mission data (Planet IDs, Asteroid ID, radii)
%
% OUTPUT:
%   kepDep[struct]      Keplerian elements of the Departure body
%   kepGA[struct]       Keplerian elements of the Gravity Assist body
%   kepAst[struct]      Keplerian elements of the Target Asteroid
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------
dateDep0 = datetime(2030,1,1,0,0,0);
mjd2000Dep = date2mjd2000(datetime2vec(dateDep0));
[kepDep] = uplanet(mjd2000Dep, data.idDep);
kepDep = kepVec2Struct(kepDep);
sDep0 = kep2car(constants.muSun,kepDep);
[~,sDep] = twoBodyProblemNoPlotHyp(constants.muSun,sDep0,1,[]);

dateGA0 = dateDep0 + 1000;
mjd2000GA = date2mjd2000(datetime2vec(dateGA0));
[kepGA] = uplanet(mjd2000GA, data.idGA);
kepGA = kepVec2Struct(kepGA);
sGA0 = kep2car(constants.muSun,kepGA);
[~,sGA] = twoBodyProblemNoPlotHyp(constants.muSun,sGA0,1,[]);

dateAst0 = datetime(2060,1,1,0,0,0);
mjd2000Ast = date2mjd2000(datetime2vec(dateAst0));
kepAst = ephAsteroids( mjd2000Ast, data.asteroidNumber);
kepAst = kepVec2Struct(kepAst);
sAst0 = kep2car(constants.muSun,kepAst);
[~,sAst] = twoBodyProblemNoPlotHyp(constants.muSun,sAst0,1,[]);

figure

scale = constants.AU;

plot3(sDep(:,1)/scale,sDep(:,2)/scale,sDep(:,3)/scale,'LineWidth',1.5)
hold on
plot3(sGA(:,1)/scale,sGA(:,2)/scale,sGA(:,3)/scale,'LineWidth',1.5)
plot3(sAst(:,1)/scale,sAst(:,2)/scale,sAst(:,3)/scale,'LineWidth',1.5)

JupCenter = [sDep(1,1)/scale,sDep(1,2)/scale,sDep(1,3)/scale];
plotPlanet('JupiterTexture.jpg',data.r_Jup*900/scale,JupCenter)

MarsCenter = [sGA(1,1)/scale,sGA(1,2)/scale,sGA(1,3)/scale];
plotPlanet('MarsTexture.jpg',data.r_Mars*8000/scale,MarsCenter)

plotPlanet('SunTexture.jpg',data.r_Sun*50/scale)


xlabel('X [AU]')
ylabel('Y [AU]')
zlabel('Z [AU]')
grid on
axis equal
legend('Dep Orbit','GA Orbit','Asteroid','Jupiter','Mars','Sun','Location','best')