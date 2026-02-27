%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORBIT PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[kepDep,kepGA,kepAst] = plotAllorbits(constants,data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESONANCE TIME ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
date_start = datetime(2030,1,1);
date_end   = datetime(2060,1,1);
date_vec   = date_start : 5 : date_end;

phase_leg1 = zeros(size(date_vec)); 
phase_leg2 = zeros(size(date_vec)); 

for k = 1:length(date_vec)
    mjd = date2mjd2000(datetime2vec(date_vec(k)));
    
    [kepDepVec] = uplanet(mjd, data.idDep);
    kepDep = kepVec2Struct(kepDepVec);
    sDep = kep2car(constants.muSun,kepDep);
    rDep = sDep(1:3);
    vDep = sDep(4:6);

   
    [kepGAVec] = uplanet(mjd, data.idGA);
    kepGA = kepVec2Struct(kepGAVec);
    sGA = kep2car(constants.muSun,kepGA);
    rGA = sGA(1:3);
    vGA = sGA(4:6);

    kepFinalArr = ephAsteroids(mjd, data.asteroidNumber);
    kepArr = kepVec2Struct(kepFinalArr);
    sArr = kep2car(constants.muSun,kepArr);  
    rAst = sArr(1:3);
    vAst = sArr(4:6);

    theta_dep = atan2(rDep(2), rDep(1));
    theta_ga = atan2(rGA(2), rGA(1));
    theta_arr = atan2(rAst(2), rAst(1));

    d_theta_leg1 = theta_ga - theta_dep;
    d_theta_leg2 = theta_arr - theta_ga;

    phase_leg1(k) = mod(rad2deg(d_theta_leg1), 360);
    phase_leg2(k) = mod(rad2deg(d_theta_leg2), 360);
end

figure('Name', 'Analisi Angoli di Fase')

subplot(2,1,1)
plot(date_vec, phase_leg1, 'b-', 'LineWidth', 1.5)
hold on 
grid on
title('Geometry Repetition: Leg 1 (Jupiter to Mars)')
axis tight
ylabel('Jupiter - Mars Phase [deg]')
ylim([0 360])
datetick('x', 'yyyy')

subplot(2,1,2)
plot(date_vec, phase_leg2, 'r-', 'LineWidth', 1.5)
hold on
grid on
title('Geometry Repetition: Leg 2 (Mars to Asteroid)')
axis tight
ylabel('Asteroid - Mars Phase [deg]')
ylim([0 360])
datetick('x', 'yyyy')

figure('Name','Pattern Recognition')
plot(date_vec, phase_leg1, 'b-', 'LineWidth', 1.3)
hold on
plot(date_vec, phase_leg2, 'r-', 'LineWidth', 1.3)
win1_start = datetime(2033,5,1);
win1_end   = datetime(2045,1,1);
win2_start = datetime(2045,1,2);
win2_end   = datetime(2056,3,1);
xr1 = xregion(win1_start, win1_end, 'FaceColor', 'g', 'FaceAlpha', 0.15, 'DisplayName', 'Ciclo 1');
xr2 = xregion(win2_start, win2_end, 'FaceColor', 'c', 'FaceAlpha', 0.15, 'DisplayName', 'Ciclo 2');
grid on
% title('Pattern Recognition ')
axis tight
ylabel('Phase [deg]')
ylim([0 360])
datetick('x', 'yyyy')
legend('Jupiter - Mars','Asteroid - Mars')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% recall synodic periods
[Tsyn1,Tsyn1_days,Tsyn1_years] = computeTsynLeg1(data.R_Jup, ...
                                                data.R_Mars,constants.muSun);

[Tsyn2,Tsyn2_days,Tsyn2_years] = computeTsynLeg2(data.R_Mars, ...
                                                kepAst.a,constants.muSun);


