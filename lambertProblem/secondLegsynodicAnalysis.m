function [vPlanet_vec,vInf_plus_vec,bestDepDate,bestArrDate,Dv_arr_mat] = secondLegsynodicAnalysis(constants, data, dateSpaceDep, dateSpaceArr, vInf_ga_max, vInf_arr_max)
% DESCRIPTION:
%   Performs the trajectory analysis for the first mission leg.
%   This function is a specialized implementation of the 'designLambertTransfer' 
%   logic.
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------

DV = NaN(length(dateSpaceDep),length(dateSpaceArr));
DV_full = NaN(length(dateSpaceDep),length(dateSpaceArr));
vInf_plus_Matrix = NaN(length(dateSpaceDep),length(dateSpaceArr),3);
vPlanet_Matrix = NaN(length(dateSpaceDep),length(dateSpaceArr),3);
ToF = NaN(length(dateSpaceDep),length(dateSpaceArr));
ToF_days = NaN(length(dateSpaceDep),length(dateSpaceArr));
sInitPlot = NaN(6,length(dateSpaceDep));
sFinalPlot = NaN(6,length(dateSpaceArr));
Dv_dep_mat = NaN(length(dateSpaceDep),length(dateSpaceArr));
Dv_arr_mat = NaN(length(dateSpaceDep),length(dateSpaceArr));

for i = 1 : length(dateSpaceDep)

    mjd2000Init = date2mjd2000(datetime2vec(dateSpaceDep(i)));

    [kepInitVec] = uplanet(mjd2000Init, data.idDep);
    kepInit = kepVec2Struct(kepInitVec);

    sInit = kep2car(constants.muSun,kepInit);

    rInit = sInit(1:3);
    vInit = sInit(4:6);
    
    for j = 1 : length(dateSpaceArr)

        mjd2000Ast = date2mjd2000(datetime2vec(dateSpaceArr(j)));
        
        kepFinalVec = ephAsteroids( mjd2000Ast, data.asteroidNumber);
        kepFinal = kepVec2Struct(kepFinalVec);

        sFinal = kep2car(constants.muSun,kepFinal);

        rFinal = sFinal(1:3);
        vFinal = sFinal(4:6);

        ToF(i,j) = seconds(dateSpaceArr(j) - dateSpaceDep(i));

        if ToF(i,j) < 0
            ToF(i,j) = NaN;
        end

        ToF_days(i,j) = ToF(i,j) / 86400;

        [~,~,~,~,vInitTras,vFinalTras,~,~] = lambertMR( rInit, rFinal, ToF(i,j), constants.muSun, 0, 0, 0 );
        
        Dv1 = norm(vInitTras' - vInit);
        Dv2 = norm(vFinal - vFinalTras');

        DV_full(i,j) = abs(Dv1) + abs(Dv2);
        Dv_dep_mat(i,j) = abs(Dv1);
        Dv_arr_mat(i,j) = abs(Dv2);

        % include constraint on Dv1 due to vInf_dep_max
        if Dv1 > vInf_ga_max
            Dv1 = NaN;
            Dv2 = NaN;
        end

        % include constraint on Dv1 due to vInf_ga_max
        if Dv2 > vInf_arr_max
            Dv1 = NaN;
            Dv2 = NaN;
        end

        DV(i,j) = abs(Dv1) + abs(Dv2);

        % load for GA analysis
        vPlanet_Matrix(i,j,:) = vInit;
        vInf_plus_Matrix(i,j,:) = vInitTras' - vInit;

        sFinalPlot(:,j) = sFinal;
    end

    sInitPlot(:,i) = sInit;
end

[DVmin,minIdx] = min(DV(:));
if isnan(DVmin)
    disp('No possile transfer in this time window with given vInf limitations')
end
[iMin, jMin] = ind2sub(size(DV), minIdx);
bestDepDate = dateSpaceDep(iMin);
bestArrDate = dateSpaceArr(jMin);
vPlanet_vec = squeeze(vPlanet_Matrix(iMin, jMin, :));
vInf_plus_vec = squeeze(vInf_plus_Matrix(iMin, jMin, :));

%% find minimum DeltaV
fprintf('\n---Second Leg Preliminary Synodic Analysis--- \n')
fprintf('minimum DV = %.4f km/s\n', DVmin);
fprintf('Best Transfer - Gravity Assist Date: %s\n', datetime(bestDepDate));
fprintf('Best Transfer - Arrival Date: %s\n', datetime(bestArrDate));

%% porkchop DeltaV plot
figure%('Color', 'w', 'Position', [100 100 1000 700])
[Xdep,Yga] = meshgrid(datenum(dateSpaceDep), datenum(dateSpaceArr));
max_Dv_visible = vInf_ga_max + vInf_arr_max;%max(vInf_dep_max,vInf_ga_max);
levels = linspace(min(DV(:)), max_Dv_visible, 40); 
contour(Xdep,Yga,DV_full',levels,'ShowText',false,"LabelFormat","%0.1f",...
    'DisplayName','\Deltav');
hold on

cb = colorbar;
cb.Label.String = '\Delta V Total  [km/s]';
% colormap(jet); 
clim([min(DV(:)), max_Dv_visible]);

xPtMin = datenum(bestDepDate);
yPtMin = datenum(bestArrDate);
plot(xPtMin,yPtMin, 'r*', 'MarkerSize', 6, 'LineWidth', 2,'DisplayName','\Deltav Min');
text(xPtMin,yPtMin, sprintf('  %.2f\n', DVmin), 'Color','r', 'FontSize',10);

contour(Xdep, Yga, Dv_dep_mat', [vInf_ga_max vInf_ga_max], ...
    'm', 'LineWidth', 1, 'DisplayName', 'GA Limit');
contour(Xdep, Yga, Dv_arr_mat', [vInf_arr_max vInf_arr_max], ...
    'r', 'LineWidth', 1, 'DisplayName', 'Arr Limit');

datetick('x','dd-mmm-yyyy','keepticks');
datetick('y','dd-mmm-yyyy','keepticks');
% xL = xlim;
% yL = ylim;
% xlim(xL);
% ylim(yL);
xlabel('Departure Date');
ylabel('Arrival Date');
xtickangle(45);
ytickangle(45);
title('Porkchop Plot: Leg 2 ');
legend('show', 'Location', 'best');
grid on

ToF_days = ToF / 86400;
ToF_years = ToF_days / 365.25;
levels_tof = floor(min(ToF_years(:))) : 0.4 : ceil(max(ToF_years(:)));
[C_tof, h_tof] = contour(Xdep, Yga, ToF_years', levels_tof, 'k','ShowText', 'on',...
    'LineWidth', 1,'DisplayName','ToF');
clabel(C_tof, h_tof, 'FontSize', 8, 'LabelSpacing', 1000,'Color', 'k');
% plot(datenum(bestDepDate), datenum(bestGADate), 'rp', 'MarkerSize', 12, 'MarkerFaceColor','w');


figure%('Color', 'w', 'Position', [100 100 1000 700]);
 
DV_clean = DV'; % Copia la matrice
mask_bad = (Dv_dep_mat > vInf_ga_max) | (Dv_arr_mat > vInf_arr_max);
DV_Clean(mask_bad) = NaN;
s = surf(Xdep, Yga, DV_clean,'HandleVisibility','off');
s.EdgeColor = 'none'; 

colormap(turbo); 
cb = colorbar;
cb.Label.String = '\Delta V Total [km/s]';

xL = xlim;
yL = ylim;
datetick('x','dd-mmm-yyyy','keepticks');
datetick('y','dd-mmm-yyyy','keepticks');
% xlim(xL);
% ylim(yL);
xlabel('Departure Date');
ylabel('Arrival Date');
xtickangle(45);
ytickangle(45);
title('Scatter Plot: Leg 2 ');
% legend('show', 'Location', 'best');
grid on


% n_total = numel(DV);
% n_bad_dep = sum(Dv_dep_mat(:) > vInf_dep_max);
% n_bad_ga  = sum(Dv_ga_mat(:) > vInf_ga_max);
% 
% fprintf('Totale Punti: %d\n', n_total);
% fprintf('Punti con Partenza Eccessiva (Bad C3): %d (%.1f%%)\n', n_bad_dep, (n_bad_dep/n_total)*100);
% fprintf('Punti con Arrivo Eccessivo (Bad GA):   %d (%.1f%%)\n', n_bad_ga, (n_bad_ga/n_total)*100);
% 
% if (n_bad_dep + n_bad_ga) == 0
%     warning('ATTENZIONE: Nessun punto viene tagliato! Controlla i valori di vInf_dep_max.');
% end

%% plot the transfer
scale = constants.AU;
ToFminDV = ToF(iMin,jMin);
time = linspace(0,ToFminDV,1e5);

% departure orbit
sInit0 = sInitPlot(:,iMin);
[tInit,sInitProp,~] = twoBodyProblemNoPlotHyp(constants.muSun,sInit0,1,[]);
[~,idxInitArr] = min(abs(tInit-ToFminDV));

% arrival orbit
sFinal0 = sFinalPlot(:,jMin);
[tFinal,sFinalProp,~] = twoBodyProblemNoPlotHyp(constants.muSun,sFinal0,-1,[]);
[~,idxFinalDep] = min(abs(tFinal+ToFminDV));

% transfer arc
rInitTras = sInit0(1:3);
rFinalTras = sFinal0(1:3);
[~,~,~,~,vInitTras,~,~,~] = lambertMR(rInitTras, rFinalTras, ToFminDV, constants.muSun, 0, 0, 0);
sTras0 = [rInitTras;vInitTras'];
[~,sTras,~] = twoBodyProblemNoPlotHyp(constants.muSun,sTras0,[],time);

figure('Name','Optimal Lambert Transfer')
% plot3(sInitPlot(1,:)/AU-0.1,sInitPlot(2,:)/AU-0.1,sInitPlot(3,:)/AU-0.1,'Color',[0.6 0.85 1.0 0.25], 'LineWidth',4)
% hold on
plot3(sInitProp(:,1)/scale,sInitProp(:,2)/scale,sInitProp(:,3)/scale,'--','color','[0.85, 0.33, 0.31]','LineWidth',1.5)
hold on
plot3(sInitProp((1:idxInitArr),1)/scale,sInitProp((1:idxInitArr),2)/scale,sInitProp((1:idxInitArr),3)/scale,'r-','LineWidth',1.5)

% plot3(sFinalPlot(1,:)/AU-0.1,sFinalPlot(2,:)/AU-0.1,sFinalPlot(3,:)/AU,'Color',[1.0 0.6 0.6 0.25], 'LineWidth',4)
plot3(sFinalProp(:,1)/scale,sFinalProp(:,2)/scale,sFinalProp(:,3)/scale,'--','color','[0.6, 0.6, 0.6]','LineWidth',1.5)
plot3(sFinalProp((1:idxFinalDep),1)/scale,sFinalProp((1:idxFinalDep),2)/scale,sFinalProp((1:idxFinalDep),3)/scale,'color','[0.6, 0.6, 0.6]','LineWidth',1.5)

plot3(sTras(:,1)/scale,sTras(:,2)/scale,sTras(:,3)/scale,'color','[0.47, 0.67, 0.19]','LineWidth',2)

% plot3(sInitProp(1,1)/scale,sInitProp(1,2)/scale,sInitProp(1,3)/scale,'bo','MarkerFaceColor','b','MarkerSize',6)
% plot3(sInitProp(idxInitArr,1)/scale,sInitProp(idxInitArr,2)/scale,sInitProp(idxInitArr,3)/scale,'bo','MarkerFaceColor','b','MarkerSize',6)
% plot3(sFinalProp(1,1)/scale,sFinalProp(1,2)/scale,sFinalProp(1,3)/scale,'ro','MarkerFaceColor','r','MarkerSize',6)
% plot3(sFinalProp(idxFinalDep,1)/scale,sFinalProp(idxFinalDep,2)/scale,sFinalProp(idxFinalDep,3)/scale,'ro','MarkerFaceColor','r','MarkerSize',6)

MarsCenter = [sInitProp(1,1)/scale,sInitProp(1,2)/scale,sInitProp(1,3)/scale];
plotPlanet('MarsTexture.jpg',data.r_Mars*8000/scale,MarsCenter)

AstCenter = [sFinalProp(1,1)/scale,sFinalProp(1,2)/scale,sFinalProp(1,3)/scale];
Gray_Dark = [0.25, 0.25, 0.25];

[x_u, y_u, z_u] = sphere(30);

X_S = x_u * 0.1 + AstCenter(1);
Y_S = y_u * 0.1 + AstCenter(2);
Z_S = z_u * 0.1 + AstCenter(3);

surf(X_S, Y_S, Z_S, ...
     'FaceColor', Gray_Dark, ...
     'EdgeColor', 'none', ...
     'FaceLighting', 'gouraud', ...
     'AmbientStrength', 0.5);

plotPlanet('SunTexture.jpg',data.r_Sun*50/scale)

grid on
axis equal
xlabel('X [AU]'); 
ylabel('Y [AU]'); 
zlabel('Z [AU]');
legend('Mars orbit','Mars motion during transfer', 'Asteroid orbit', ...
    'Asteroid motion during transfer', 'Transfer arc', 'Location','best')

end