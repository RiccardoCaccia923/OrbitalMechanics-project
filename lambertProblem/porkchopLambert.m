function [DV, iMin, jMin, ToF, dateSpaceInit, dateSpaceFinal, sInitPlot, sFinalPlot] = ...
            porkchopLambert(idInit, idFinal, dateInit0, dateInitEnd, dateFinal0, dateFinalEnd, dateStep, vInf)
%
% DESCRIPTION:
%   Performs a systematic grid search over specified departure and arrival 
%   windows to solve Lambert's problem. 
%   The function computes the total Delta V for each pair of dates and 
%   generates a "Porkchop Plot" (contour map of mission cost). 
%   It automatically identifies and highlights the global minimum Delta V 
%   solution within the search space.
%
% INPUT:
%   idInit[1]       ID of the departure body
%   idFinal[1]      ID of the arrival body
%   dateInit0       Start of departure window (datetime)
%   dateInitEnd     End of departure window (datetime)
%   dateFinal0      Start of arrival window (datetime)
%   dateFinalEnd    End of arrival window (datetime)
%   dateStep        Time step for the grid search (duration)
%   vInf[1]         Upper limit for Delta V contour levels in the plot [km/s]
%
% OUTPUT:
%   DV[NxM]         Matrix of total Delta V costs [km/s]
%   iMin, jMin      Indices of the minimum cost solution
%   ToF[NxM]        Time of Flight matrix [s]
%   dateSpace...    Vectors of dates used for the axes
%   s...Plot        State vectors of bodies for plotting purposes
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------

% define constants
muSun = astroConstants(4);

% temporal grid
dateSpaceInit  = dateInit0:dateStep:dateInitEnd;
dateSpaceFinal = dateFinal0:dateStep:dateFinalEnd;

DV   = zeros(length(dateSpaceInit), length(dateSpaceFinal));
ToF  = zeros(length(dateSpaceInit), length(dateSpaceFinal));
sInitPlot  = zeros(6, length(dateSpaceInit));
sFinalPlot = zeros(6, length(dateSpaceFinal));

% cycle on departure dates
for i = 1:length(dateSpaceInit)

    mjd2000Init = date2mjd2000(datetime2vec(dateSpaceInit(i)));
    kepInitVec  = uplanet(mjd2000Init, idInit);
    kepInit     = kepVec2Struct(kepInitVec);
    sInit       = kep2car(muSun, kepInit);

    rInit = sInit(1:3);
    vInit = sInit(4:6);

    % cycle on arrival dates
    for j = 1:length(dateSpaceFinal)

        mjd2000Final = date2mjd2000(datetime2vec(dateSpaceFinal(j)));
        kepFinalVec  = uplanet(mjd2000Final, idFinal);
        kepFinal     = kepVec2Struct(kepFinalVec);
        sFinal       = kep2car(muSun, kepFinal);

        rFinal = sFinal(1:3);
        vFinal = sFinal(4:6);

        ToF(i,j) = seconds(dateSpaceFinal(j) - dateSpaceInit(i));

        [~,~,~,~,vInitTras,vFinalTras,~,~] = lambertMR(rInit, rFinal, ToF(i,j), muSun, 0, 0, 0);

        Dv1 = norm(vInitTras' - vInit);
        Dv2 = norm(vFinal - vFinalTras');
        DV(i,j) = abs(Dv1) + abs(Dv2);

        sFinalPlot(:,j) = sFinal;
    end

    sInitPlot(:,i) = sInit;
end

%% minimum DV
[DVmin, minIdx] = min(DV(:));
[iMin, jMin]    = ind2sub(size(DV), minIdx);
bestDepDate     = dateSpaceInit(iMin);
bestArrDate     = dateSpaceFinal(jMin);

fprintf('minimum DV = %.4f km/s\n', DVmin);
fprintf('Best Transfer - Departure Date: %s\n', datetime(bestDepDate));
fprintf('Best Transfer - Arrival Date : %s\n', datetime(bestArrDate));

%% Porkchop Plot
[Xdep, Yarr] = meshgrid(datenum(dateSpaceInit),datenum(dateSpaceFinal));

figure
levels = 5:1:vInf;
[C,h] = contour(Xdep,Yarr,DV',levels,'ShowText',true,"LabelFormat","%0.1f");
clabel(C, h, 'LabelSpacing', 300);
colorbar
clim([5 10]);
hold on;
xPtMin = datenum(bestDepDate);
yPtMin = datenum(bestArrDate);
plot(xPtMin, yPtMin, 'r*', 'MarkerSize', 6, 'LineWidth', 2);
text(xPtMin, yPtMin, sprintf('  %.2f', DVmin), 'Color','r', 'FontSize',9);
hold off;
datetick('x','dd-mmm-yyyy','keepticks');
datetick('y','dd-mmm-yyyy','keepticks');
xlabel('Departure Date');
ylabel('Arrival Date');
xtickangle(45);
ytickangle(45);
title('\DeltaV [km/s]');

end
