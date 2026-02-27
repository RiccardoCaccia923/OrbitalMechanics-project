%% backup script for synodic analysis

dateDep0 = datetime(2030,1,1,0,0,0);
% to be conservative and accounting for the fact that in reality orbits are
% elliptic, the computed Tsyn will be augumented with a safety buffer
dateDepEnd = dateDep0 + Tsyn_days*1.2;

% since this is only a preliminary analysis it can be intersting to search
% for solutions with a ToF lower than the one computed for the hohmann
% transfer, which is known to be maybe the cheapest possible transfer but
% for sure not the shortest one but also consider for tof higher than the
% hohmann to be conservative
dateGA0 = dateDep0 + 0.2*days(hohmann.TOF_days);
dateGAEnd = dateDepEnd + 2*days(hohmann.TOF_days);

% define the time window vectors
dateNum = 200;
dateSpaceDep = linspace(dateDep0,dateDepEnd,dateNum);
dateSpaceGA = linspace(dateGA0,dateGAEnd,dateNum);

DV = NaN(length(dateSpaceDep),length(dateSpaceGA));
ToF = NaN(length(dateSpaceDep),length(dateSpaceGA));
sInitPlot = NaN(6,length(dateSpaceDep));
sFinalPlot = NaN(6,length(dateSpaceGA));

for i = 1 : length(dateSpaceDep)

    mjd2000Init = date2mjd2000(datetime2vec(dateSpaceDep(i)));

    [kepInitVec] = uplanet(mjd2000Init, data.idInit);
    kepInit = kepVec2Struct(kepInitVec);

    sInit = kep2car(constants.muSun,kepInit);

    rInit = sInit(1:3);
    vInit = sInit(4:6);
    
    for j = 1 : length(dateSpaceGA)

        mjd2000Final = date2mjd2000(datetime2vec(dateSpaceGA(j)));

        [kepFinalVec] = uplanet(mjd2000Final, data.idFinal);
        kepFinal = kepVec2Struct(kepFinalVec);

        sFinal = kep2car(constants.muSun,kepFinal);

        rFinal = sFinal(1:3);
        vFinal = sFinal(4:6);

        ToF(i,j) = seconds(dateSpaceGA(j) - dateSpaceDep(i));

        if ToF(i,j) < 0
            ToF(i,j) = NaN;
        end

        [~,~,~,~,vInitTras,vFinalTras,~,~] = lambertMR( rInit, rFinal, ToF(i,j), constants.muSun, 0, 0, 0 );
        
        Dv1 = norm(vInitTras' - vInit);
        Dv2 = norm(vFinal - vFinalTras');

        DV_full(i,j) = abs(Dv1) + abs(Dv2);
        Dv_dep_mat(i,j) = abs(Dv1);
        Dv_ga_mat(i,j) = abs(Dv2);

        % include constraint on Dv1 due to vInf_dep_max
        if Dv1 > vInf_dep_max
            Dv1 = NaN;
            Dv2 = NaN;
        end

        % include constraint on Dv1 due to vInf_ga_max
        if Dv2 > vInf_ga_max
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
bestDepDate = dateSpaceDep(iMin);
bestGADate = dateSpaceGA(jMin);


%% find minimum DeltaV
fprintf('\n---First Leg Preliminary Synodic Analysis--- \n')
fprintf('minimum DV = %.4f km/s\n', DVmin);
fprintf('Best Transfer - Departure Date: %s\n', datetime(bestDepDate));
fprintf('Best Transfer - Gravity Assist Date : %s\n', datetime(bestGADate));

%% porkchop DeltaV plot
figure('Color', 'w', 'Position', [100 100 1000 700])
[Xdep,Yga] = meshgrid(datenum(dateSpaceDep), datenum(dateSpaceGA));
max_Dv_visible = vInf_dep_max + vInf_ga_max;%max(vInf_dep_max,vInf_ga_max);
levels = linspace(min(DV(:)), max_Dv_visible, 40); 
[C,h] = contour(Xdep,Yga,DV_full',levels,'ShowText',false,"LabelFormat","%0.1f",...
    'DisplayName','\Deltav');
hold on
%clabel(C, h, 'LabelSpacing', 400)
cb = colorbar;
cb.Label.String = '\Delta V Total  [km/s]';
% colormap(jet); 
clim([min(DV(:)), max_Dv_visible]);

xPtMin = datenum(bestDepDate);
yPtMin = datenum(bestGADate);
plot(xPtMin,yPtMin, 'r*', 'MarkerSize', 6, 'LineWidth', 2,'DisplayName','\Deltav Min');
text(xPtMin,yPtMin, sprintf('  %.2f\n', DVmin), 'Color','r', 'FontSize',10);

contour(Xdep, Yga, Dv_dep_mat', [vInf_dep_max vInf_dep_max], ...
    'm', 'LineWidth', 1, 'DisplayName', 'C3 Limit');
contour(Xdep, Yga, Dv_ga_mat', [vInf_ga_max vInf_ga_max], ...
    'r', 'LineWidth', 1, 'DisplayName', 'GA Limit');

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
title('Porkchop Plot: Leg 1 (Contours = Constraints)');
legend('show', 'Location', 'best');
grid on

ToF_days = ToF / 86400;
ToF_years = ToF_days / 365.25;
levels_tof = floor(min(ToF_years(:))) : 0.4 : ceil(max(ToF_years(:)));
[C_tof, h_tof] = contour(Xdep, Yga, ToF_years', levels_tof, 'k','ShowText', 'on',...
    'LineWidth', 1,'DisplayName','ToF');
clabel(C_tof, h_tof, 'FontSize', 8, 'LabelSpacing', 1000,'Color', 'k');
% plot(datenum(bestDepDate), datenum(bestGADate), 'rp', 'MarkerSize', 12, 'MarkerFaceColor','w');


figure('Color', 'w', 'Position', [100 100 1000 700]);
 
DV_clean = DV'; % Copia la matrice
mask_bad = (Dv_dep_mat > vInf_dep_max) | (Dv_ga_mat > vInf_ga_max);
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
title('Porkchop Plot: Leg 1 (Contours = Constraints)');
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
