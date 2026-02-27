function [DV, iMin, jMin, ToF, dateSpaceInit, dateSpaceFinal, sInitPlot, sFinalPlot] = ...
            designLambertTransfer(idInit, idFinal, dateInit0, dateInitEnd, dateFinal0, dateFinalEnd, dateStep, dateNum, vInf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% designLambertTransfer  Optimizes and visualizes an interplanetary transfer.
%
% PROTOTYPE:
%   [DV, iMin, jMin, ToF, dateSpaceInit, dateSpaceFinal, sInitPlot, sFinalPlot] = ...
%       designLambertTransfer(idInit, idFinal, dateInit0, dateInitEnd, ...
%                             dateFinal0, dateFinalEnd, dateStep, dateNum, vInf)
%
% DESCRIPTION:
%   Performs a grid search over specified departure and arrival windows to 
%   find the optimal interplanetary transfer trajectory between two bodies 
%   (e.g., Earth and Mars) using the Lambert solver.
%   The function:
%   1. Discretizes the date ranges (either by step size or number of points).
%   2. Solves Lambert's Problem for each (Departure, Arrival) pair.
%   3. Computes the total Delta V (Departure + Arrival) and filters out 
%      solutions exceeding the maximum allowable hyperbolic excess velocity (vInf).
%   4. Identifies the global minimum Delta V solution.
%   5. Generates a Porkchop Plot (Delta V contours) and a 3D trajectory plot.
%
% INPUT:
%   idInit[1]       NAIF ID of the departure body (e.g., 399 for Earth)
%   idFinal[1]      NAIF ID of the arrival body (e.g., 499 for Mars)
%   dateInit0       Start of departure window (datetime or datenum)
%   dateInitEnd     End of departure window (datetime or datenum)
%   dateFinal0      Start of arrival window (datetime or datenum)
%   dateFinalEnd    End of arrival window (datetime or datenum)
%   dateStep        Time step for grid search (duration, e.g., days(1)); 
%                   leave empty [] if using dateNum
%   dateNum[1]      Number of points for grid discretization; 
%                   leave empty [] if using dateStep
%   vInf[1]         Maximum allowable hyperbolic excess velocity [km/s]
%                   (Used to filter feasible departure options)
%
% OUTPUT:
%   DV[NxM]         Matrix of total Delta V costs [km/s]
%   iMin[1]         Row index of the minimum Delta V solution
%   jMin[1]         Column index of the minimum Delta V solution
%   ToF[NxM]        Matrix of Time of Flight values [s]
%   dateSpaceInit   Vector of departure dates used in the grid
%   dateSpaceFinal  Vector of arrival dates used in the grid
%   sInitPlot       State vectors of the departure body over the window
%   sFinalPlot      State vectors of the arrival body over the window
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------
muSun = astroConstants(4);

if isempty(dateStep) && ~isempty(dateNum)
    dateSpaceInit = linspace(dateInit0,dateInitEnd,dateNum);
    dateSpaceFinal = linspace(dateFinal0,dateFinalEnd,dateNum);
elseif isempty(dateNum) && ~isempty(dateStep)
    dateSpaceInit = dateInit0:dateStep:dateInitEnd;
    dateSpaceFinal = dateFinal0:dateStep:dateFinalEnd;
else
    disp('cant define date vector')
end

DV = zeros(length(dateSpaceInit),length(dateSpaceFinal));
ToF = zeros(length(dateSpaceInit),length(dateSpaceFinal));
sInitPlot = zeros(6,length(dateSpaceInit));
sFinalPlot = zeros(6,length(dateSpaceFinal));

for i = 1 : length(dateSpaceInit)

    mjd2000Init = date2mjd2000(datetime2vec(dateSpaceInit(i)));

    [kepInitVec] = uplanet(mjd2000Init, idInit);
    kepInit = kepVec2Struct(kepInitVec);

    sInit = kep2car(muSun,kepInit);

    rInit = sInit(1:3);
    vInit = sInit(4:6);
    
    for j = 1 : length(dateSpaceFinal)

        mjd2000Final = date2mjd2000(datetime2vec(dateSpaceFinal(j)));

        [kepFinalVec] = uplanet(mjd2000Final, idFinal);
        kepFinal = kepVec2Struct(kepFinalVec);

        sFinal = kep2car(muSun,kepFinal);

        rFinal = sFinal(1:3);
        vFinal = sFinal(4:6);

        ToF(i,j) = seconds(dateSpaceFinal(j) - dateSpaceInit(i));

        [~,~,~,~,vInitTras,vFinalTras,~,~] = lambertMR( rInit, rFinal, ToF(i,j), muSun, 0, 0, 0 );
        
        Dv1 = norm(vInitTras' - vInit);
        Dv2 = norm(vFinal - vFinalTras');

        % include constraint on Dv1 due to Vinf
        if Dv1 > vInf
            Dv1 = NaN;
            Dv2 = NaN;
        end

        DV(i,j) = abs(Dv1) + abs(Dv2);

        sFinalPlot(:,j) = sFinal;
    end

    sInitPlot(:,i) = sInit;
end

[DVmin,minIdx] = min(DV(:));
if isnan(DVmin)
    disp('No possile transfer in this time window with given Vinf')
end
[iMin, jMin] = ind2sub(size(DV), minIdx);
bestDepDate = dateSpaceInit(iMin);
bestArrDate = dateSpaceFinal(jMin);

%% find minimum DeltaV
fprintf('minimum DV = %.4f km/s\n', DVmin);
fprintf('Best Transfer - Departure Date: %s\n', datetime(bestDepDate));
fprintf('Best Transfer - Arrival Date : %s\n', datetime(bestArrDate));
[Xdep,Yarr] = meshgrid(datenum(dateSpaceInit), datenum(dateSpaceFinal));

%% porkchop DeltaV plot
figure
% levels = 5:1:vInf;   
[C,h] = contour(Xdep,Yarr,DV','ShowText',true,"LabelFormat","%0.1f");
clabel(C, h, 'LabelSpacing', 400);
colorbar
% clim([5 10]);
hold on;
xPtMin = datenum(dateSpaceInit(iMin));
yPtMin = datenum(dateSpaceFinal(jMin));
plot(xPtMin, yPtMin, 'r*', 'MarkerSize', 6, 'LineWidth', 2);
text(xPtMin, yPtMin, sprintf('  %.2f', DVmin), 'Color','r', 'FontSize',9);

% ToFdays = ToF/86400;
% ToFvalues = [60 120 180 240 300];
% [CTF, hTF] = contour(Xdep, Yarr, ToFdays', ToFvalues,'LineColor','k','LineWidth',1);
% clabel(CTF, hTF, 'FontWeight','bold');
% hold off;
xL = xlim;
yL = ylim;
datetick('x','dd-mmm-yyyy','keepticks');
datetick('y','dd-mmm-yyyy','keepticks');
xlim(xL);
ylim(yL);
xlabel('Departure Date');
ylabel('Arrival Date');
xtickangle(45);
ytickangle(45);
title('\DeltaV [km/s]');

%% refine the solution
%fminunc()

%% plot the transfer
AU  = 149597870.7;
ToFminDV = ToF(iMin,jMin);
time = linspace(0,ToFminDV,1e5);

% departure orbit
sInit0 = sInitPlot(:,iMin);
[tInit,sInitProp,~] = twoBodyProblemNoPlot(muSun,sInit0,1,[]);
[~,idxInitArr] = min(abs(tInit-ToFminDV));

% arrival orbit
sFinal0 = sFinalPlot(:,jMin);
[tFinal,sFinalProp,~] = twoBodyProblemNoPlot(muSun,sFinal0,-1,[]);
[~,idxFinalDep] = min(abs(tFinal+ToFminDV));

% transfer arc
rInitTras = sInit0(1:3);
rFinalTras = sFinal0(1:3);
[~,~,~,~,vInitTras,vFinalTras,~,~] = lambertMR(rInitTras, rFinalTras, ToFminDV, muSun, 0, 0, 0);
sTras0 = [rInitTras;vInitTras'];
[tTras,sTras,~] = twoBodyProblemNoPlot(muSun,sTras0,[],time);

figure('Name','Interplanetary Lambert Transfer')
% plot3(sInitPlot(1,:)/AU-0.1,sInitPlot(2,:)/AU-0.1,sInitPlot(3,:)/AU-0.1,'Color',[0.6 0.85 1.0 0.25], 'LineWidth',4)
% hold on
plot3(sInitProp(:,1)/AU,sInitProp(:,2)/AU,sInitProp(:,3)/AU,'b--','LineWidth',2)
hold on
plot3(sInitProp((1:idxInitArr),1)/AU,sInitProp((1:idxInitArr),2)/AU,sInitProp((1:idxInitArr),3)/AU,'b-','LineWidth',2)

% plot3(sFinalPlot(1,:)/AU-0.1,sFinalPlot(2,:)/AU-0.1,sFinalPlot(3,:)/AU,'Color',[1.0 0.6 0.6 0.25], 'LineWidth',4)
plot3(sFinalProp(:,1)/AU,sFinalProp(:,2)/AU,sFinalProp(:,3)/AU,'r--','LineWidth',2)
plot3(sFinalProp((1:idxFinalDep),1)/AU,sFinalProp((1:idxFinalDep),2)/AU,sFinalProp((1:idxFinalDep),3)/AU,'r-','LineWidth',2)

plot3(sTras(:,1)/AU,sTras(:,2)/AU,sTras(:,3)/AU,'g-','LineWidth',2)

plot3(sInitProp(1,1)/AU,sInitProp(1,2)/AU,sInitProp(1,3)/AU,'bo','MarkerFaceColor','b','MarkerSize',6)
plot3(sInitProp(idxInitArr,1)/AU,sInitProp(idxInitArr,2)/AU,sInitProp(idxInitArr,3)/AU,'bo','MarkerFaceColor','b','MarkerSize',6)
plot3(sFinalProp(1,1)/AU,sFinalProp(1,2)/AU,sFinalProp(1,3)/AU,'ro','MarkerFaceColor','r','MarkerSize',6)
plot3(sFinalProp(idxFinalDep,1)/AU,sFinalProp(idxFinalDep,2)/AU,sFinalProp(idxFinalDep,3)/AU,'ro','MarkerFaceColor','r','MarkerSize',6)

plot3(0,0,0,'yo','MarkerFaceColor','y','MarkerSize',10)

grid on
axis equal
xlabel('x [AU]'); 
ylabel('y [AU]'); 
zlabel('z [AU]');
legend('Earth orbit','Earth motion during transfer', 'Planet orbit', ...
    'Planet motion during transfer', 'Transfer arc', 'Location','best')


end