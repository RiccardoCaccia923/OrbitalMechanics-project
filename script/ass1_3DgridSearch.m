%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- 3D GRID SEARCH --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define dates

% first window
% dateDep0 = datetime(2031,5,1);
% dateDepEnd = datetime(2040,1,1);
% 
% dateGA0 = datetime(2035,1,1);
% dateGAEnd = datetime(2043,1,1);
% 
% dateArr0 = datetime(2036,1,1);
% dateArrEnd = datetime(2045,1,1);

% second window
% dateDep0 = datetime(2043,1,1);
% dateDepEnd = datetime(2050,1,1);
% 
% dateGA0 = datetime(2047,1,1);
% dateGAEnd = datetime(2053,1,1);
% 
% dateArr0 = datetime(2048,1,1);
% dateArrEnd = datetime(2056,1,1);

% narrow window
dateDep0 = datetime(2044,1,1);
dateDepEnd = datetime(2045,1,1);

dateGA0 = datetime(2047,1,1);
dateGAEnd = datetime(2048,1,1);

dateArr0 = datetime(2049,1,1);
dateArrEnd = datetime(2050,1,1);


dateNum = 100;
dateSpaceDep = linspace(dateDep0,dateDepEnd,dateNum);
dateSpaceGA = linspace(dateGA0,dateGAEnd,dateNum);
dateSpaceArr = linspace(dateArr0,dateArrEnd,dateNum);

% pre allocation
Dv_dep = NaN(length(dateSpaceDep),length(dateSpaceGA));
Dv_arr = NaN(length(dateSpaceGA),length(dateSpaceArr));
Dv_ga = NaN(length(dateSpaceDep), length(dateSpaceGA), length(dateSpaceArr));

tof1 = NaN(length(dateSpaceDep),length(dateSpaceGA));
tof1_days = NaN(length(dateSpaceDep),length(dateSpaceGA));
tof2 = NaN(length(dateSpaceGA),length(dateSpaceArr));
tof2_days = NaN(length(dateSpaceGA),length(dateSpaceArr));

vInf_in_mat = NaN(length(dateSpaceDep),length(dateSpaceGA),3);
vInf_out_mat = NaN(length(dateSpaceGA),length(dateSpaceArr),3);

J_total = zeros(length(dateSpaceDep), length(dateSpaceGA), length(dateSpaceArr));

% 3 nested loops
for i = 1:length(dateSpaceDep)
    mjd2000Dep = date2mjd2000(datetime2vec(dateSpaceDep(i)));
    [kepDepVec] = uplanet(mjd2000Dep, data.idDep);
    kepDep = kepVec2Struct(kepDepVec);
    sDep = kep2car(constants.muSun,kepDep);
    rDep = sDep(1:3);
    vDep = sDep(4:6);

    for j = 1:length(dateSpaceGA)
        mjd2000GA = date2mjd2000(datetime2vec(dateSpaceGA(j)));
        [kepGAVec] = uplanet(mjd2000GA, data.idGA);
        kepGA = kepVec2Struct(kepGAVec);
        sGA = kep2car(constants.muSun,kepGA);
        rGA = sGA(1:3);
        vGA = sGA(4:6);

        tof1(i,j) = seconds(dateSpaceGA(j) - dateSpaceDep(i));
        tof1_days(i,j) = tof1(i,j)/86400;

        [~,~,~,~, v1_dep, v1_arr, ~, ~] = lambertMR(rDep, rGA, tof1(i,j), constants.muSun, 0, 0, 0);
        
        Dv_dep(i,j) = norm(v1_dep' - vDep);
        vInf_in = v1_arr' - vGA;
        vInf_in_mat(i,j,:) = vInf_in;


        for k = 1:length(dateSpaceArr)
            mjd2000Ast = date2mjd2000(datetime2vec(dateSpaceArr(k)));
            kepFinalArr = ephAsteroids( mjd2000Ast, data.asteroidNumber);
            kepArr = kepVec2Struct(kepFinalArr);
            sArr = kep2car(constants.muSun,kepArr);  
            rAst = sArr(1:3);
            vAst = sArr(4:6);

            tof2(j,k) = seconds(dateSpaceArr(k) - dateSpaceGA(j));
            tof2_days(j,k) = tof2(j,k)/86400;

            [~,~,~,~, v2_dep, v2_arr, ~, ~] = lambertMR(rGA, rAst, tof2(j,k), constants.muSun, 0, 0, 0);
             
            vInf_out = v2_dep' - vGA;
            vInf_out_mat(j,k,:) = vInf_out;              
            % Dv_ga(i,j,k) = norm(vInf_out - vInf_in);
            v_in_mag = norm(vInf_in);
            v_out_mag = norm(vInf_out);
            cos_delta = dot(vInf_in, vInf_out) / (v_in_mag * v_out_mag);
            delta_req = acos(max(-1, min(1, cos_delta)));
            
            % 3. Angolo massimo fornito dalla gravità a r_safe
            r_safe = data.r_Mars + data.h_atmMars;    
            e_in  = 1 + (r_safe * v_in_mag^2) / data.muMars;
            e_out = 1 + (r_safe * v_out_mag^2) / data.muMars;
            delta_grav = asin(1/e_in) + asin(1/e_out);

            if delta_req <= delta_grav
                Dv_ga(i,j,k) = abs(v_out_mag - v_in_mag);
            else

                delta_maneuver = delta_req - delta_grav;
                vp_in  = sqrt(v_in_mag^2 + 2*data.muMars/r_safe);
                vp_out = sqrt(v_out_mag^2 + 2*data.muMars/r_safe);
                Dv_ga(i,j,k) = sqrt(vp_in^2 + vp_out^2 - 2*vp_in*vp_out*cos(delta_maneuver));
            end 

            Dv_arr(j,k) = norm(vAst - v2_arr');     
             
            % total mission cost
            J_total(i,j,k) = Dv_dep(i,j) + Dv_ga(i,j,k) + Dv_arr(j,k);
            J_full(i,j,k) = Dv_dep(i,j) + Dv_ga(i,j,k) + Dv_arr(j,k);

            % constraints on departure and arrival
            if Dv_dep(i,j) > Leg1.vInf_dep_max || Dv_arr(j,k) > Leg2.vInf_arr_max
                J_total(i,j,k) = NaN;
            end

            % constraint on gravity assist
            % delta_required = acos(dot(vInf_in, vInf_out) / (norm(vInf_in)*norm(vInf_out)));
            % rp_safe = data.r_Mars + data.h_atmMars;
            % e_in_safe  = 1 + (rp_safe * norm(vInf_in)^2) / data.muMars;
            % e_out_safe = 1 + (rp_safe * norm(vInf_out)^2) / data.muMars;
            % delta_max_possible = asin(1/e_in_safe) + asin(1/e_out_safe);
            % if delta_required > delta_max_possible
            % J_total(i,j,k) = NaN;
            % end
        end
    end
end

% find optimal transfer
[min_cost, linear_idx] = min(J_total(:));

[i_opt, j_opt, k_opt] = ind2sub(size(J_total), linear_idx);

best_DepDate = dateSpaceDep(i_opt);
best_GADate  = dateSpaceGA(j_opt);
best_ArrDate = dateSpaceArr(k_opt);

opt_Dv_dep = Dv_dep(i_opt, j_opt);
opt_Dv_ga  = Dv_ga(i_opt, j_opt, k_opt);
opt_Dv_arr = Dv_arr(j_opt, k_opt);

fprintf('\n---OPTIMAL SOLUTION (GRID SEARCH)---\n')
fprintf('Jupiter Departure:  %s\n', datetime(best_DepDate))
fprintf('Mars Flyby:         %s\n', datetime(best_GADate))
fprintf('Asteroid Arrival:   %s\n', datetime(best_ArrDate))
fprintf('--------------------------------------\n')
fprintf('Total Delta V:      %.4f km/s\n', min_cost)
fprintf('  - Departure:      %.4f km/s\n', opt_Dv_dep)
fprintf('  - Flyby (GA):     %.4f km/s\n', opt_Dv_ga)
fprintf('  - Arrival:        %.4f km/s\n', opt_Dv_arr)
fprintf('--------------------------------------\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT COMPLETE TRANSFER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot3Opt(constants,data,best_DepDate,best_GADate,best_ArrDate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAVITY ASSIST ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt_vInf_in = squeeze(vInf_in_mat(i_opt,j_opt,:));
opt_vInf_out = squeeze(vInf_out_mat(j_opt,k_opt,:));
GA = computeGA(data,opt_vInf_in,opt_vInf_out);
plotGAOpt(data,GA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PORKCHOP PLOT ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_vec = datenum(dateSpaceDep);
y_vec = datenum(dateSpaceArr); 

% LOCAL PORKCHOP
Z_slice = squeeze(J_full(:, j_opt, :)); 

figure('Name', 'Porkchop Plot: Fixed Flyby', 'Color', 'w')
[C, h] = contourf(x_vec, y_vec, Z_slice', 50, 'LineColor', 'none');
hold on
grid on

% colormap(jet) 
c = colorbar; 
c.Label.String = 'Delta V Totale [km/s]';
clim([min(J_total(:)), min(J_total(:)) + 2]); % Mostra solo i 5 km/s sopra il minimo

plot(x_vec(i_opt), y_vec(k_opt), 'p', 'MarkerFaceColor', 'w', 'MarkerSize', 14)

datetick('x', 'mmm-yy', 'keeplimits')
datetick('y', 'mmm-yy', 'keeplimits')
xlabel('Data Partenza (Giove)')
ylabel('Data Arrivo (Asteroide)')
% title(['Porkchop Plot (Flyby fissato al ' datenum(dateSpaceGA(j_opt)) ')'])


% GLOBAL PORCKCHOP
% min(J_total, [], 2) 
[Z_global, best_GA_indices] = min(J_full, [], 2); 
Z_global = squeeze(Z_global); % Diventa 2D (nDep x nArr)

figure('Name', 'Porkchop Plot: Optimized Flyby', 'Color', 'w');
[C2, h2] = contourf(x_vec, y_vec, Z_global', 50, 'LineColor', 'none');
hold on; grid on;

% colormap(jet); 
c = colorbar; 
c.Label.String = 'Total Delta V (with best GA) [km/s]';
clim([min(J_total(:)), min(J_total(:)) + 10]); % Focus sulle zone buone

plot(x_vec(i_opt), y_vec(k_opt), 'p', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 15);

datetick('x', 'mmm-yy', 'keeplimits');
datetick('y', 'mmm-yy', 'keeplimits');
xtickangle(45);
xlabel('Departure Date (Jupiter)');
ylabel('Arrival Date (Asteroid)');
% title('Global Mission Porkchop');

[X_mesh, Y_mesh] = meshgrid(x_vec, y_vec);
ToF_years = (Y_mesh - X_mesh)/365.25; 
[C_tof, h_tof] = contour(X_mesh, Y_mesh, ToF_years, 10, 'k', 'ShowText', 'on');
h_tof.LevelList = round(h_tof.LevelList); 
clabel(C_tof, h_tof, 'Color', 'k', 'FontSize', 8);