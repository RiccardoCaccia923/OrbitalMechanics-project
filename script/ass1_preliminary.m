%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESIGN of GRAVITY-ASSISTED INTERPLANETARY TRANSFER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONSTANTS
constants.muSun = astroConstants(4);
constants.AU = astroConstants(2);

%% DATA
data.r_Sun = 696340;
data.idDep = 5;
data.idGA = 4;
data.asteroidNumber = 448628;
data.R_Jup  = 5.204 * constants.AU;  
data.r_Jup = 71492;
data.R_Mars = 1.524 * constants.AU;  
data.muMars = astroConstants(14);
data.r_Mars = 3396;
data.h_atmMars = 120;

%% ORBIT PLOT
[kepDep,kepGA,kepAst] = plotAllorbits(constants,data);

%% PRELIMINARY ANALYSIS
% in this problem there are 3 DoF, so 3D grid search is not an intelligent
% approach due to the fact that it is very computationally expensive.
% the idea is to start by trying to implement an approach based on two 2d
% grid search conbined.
% an important aspect is to reduce the dimentions of time windows and this 
% will be done by exploiting different considerations

% 1. hohmann transfer between two orbit usually is the cheapest transfer,
%    so it can show the minimum value of Dv for the transfer and also gives 
%    the tof which minimizes the Dv so i will find:
%    - near minimun DeltaV
%    - an high but efficient ToF 
[Leg1.hohmann.dV, Leg1.hohmann.TOF, Leg1.hohmann.TOF_days] = ...
    hohmannPlanets(data.R_Jup,data.R_Mars);
                                                    
% 2. synodic period rapresents the time which is passed from last time in
%    which two bodies were in the same relative position, so by analyzing
%    a full synodic period it is possible to evaluate all possible
%    configurations (obviously only for the 2 body involed)
%    in this way computing over 1 synodic period the minimum tof for the
%    two different tbp i can determine limit fot the time window of tGA
%    
%    before doing Tsyn analysis a crucial aspect is to define limit values
%    for DeltaV, in order to clean the grid from non feasible cases
%    DeltaV limits are different if we consider the departure or the
%    arrival planer, so:
%    - departure: launcer limit
%    - arrival: maximum turning angle*
%
%    * since at the arrival planet a gravity assist manoeuve must be
%      performed, is important to understand that if the velocity of
%      arrival is too high respect to the gravity acceleration of that
%      planet the manoeuvre won't happen
%    
%    both this values must be conservative, to not compromise the whole
%    analysis, so a safe margin factor will be applied to both

% define first leg constraints
 constraints.rocketC3 = 50;                                                
 constraints.deltaMin = deg2rad(5);                                                

% compute vInf limitations from leg constraints
[Leg1.vInf_dep_max,Leg1.vInf_ga_max] = vInfLimitLeg1(data,constraints);

% before starting iterations lets define synodic period
[Leg1.Tsyn,Leg1.Tsyn_days,Leg1.Tsyn_years] = computeTsynLeg1(data.R_Jup, ...
                                                data.R_Mars,constants.muSun);


%% First Leg Preliminary analysis
dateDep0 = datetime(2044,1,1,0,0,0);
% to be conservative and accounting for the fact that in reality orbits are
% elliptic, the computed Tsyn will be augumented with a safety buffer
dateDepEnd = dateDep0 + Leg1.Tsyn_days*1.2;
% since this is only a preliminary analysis it can be intersting to search
% for solutions with a ToF lower than the one computed for the hohmann
% transfer, which is known to be maybe the cheapest possible transfer but
% for sure not the shortest one but also consider for tof higher than the
% hohmann to be conservative
% dateGA0 = datetime(2035,8,15,0,0,0)-100;
% dateGAEnd = datetime(2035,5,15,0,0,0)+1000;
dateGA0 = dateDep0 + 0.5*days(Leg1.hohmann.TOF_days);
dateGAEnd = dateDepEnd + 5*days(Leg1.hohmann.TOF_days);
% define the time window dimension
dateNum = 100;
dateSpaceDep = linspace(dateDep0,dateDepEnd,dateNum);
dateSpaceGA1 = linspace(dateGA0,dateGAEnd,dateNum);

[GA.V_minus,GA.vInf_minus,...
    Leg1.BestDepDate,Leg1.bestGAdate,GA.matDvDep] = firstLegsynodicAnalysis(constants, data,...
                                                        dateSpaceDep, dateSpaceGA1, ...
                                                        Leg1.vInf_dep_max, Leg1.vInf_ga_max);

%% Second Leg Preliminary Analysis
% the idea is to repeat the same analysis done for the first leg wit the
% second

% estimate with hohmann only the ToF, so even if a change of plan would be
% necessary it won't be considered since it is an instantaneous manoeuvre
% so it doesn't change the ToF
[Leg2.hohmann.TOF, Leg2.hohmann.TOF_days] = hohmannAsteroidToF(data.R_Mars, ...
                                                kepAst.a);

% compute the synodic period between the two bodies
[Leg2.Tsyn,Leg2.Tsyn_days,Leg2.Tsyn_years] = computeTsynLeg2(data.R_Mars,kepAst.a,constants.muSun);

% define vInf limitations
[Leg2.vInf_ga_max,Leg2.vInf_arr_max] = vInfLimitLeg2(data,constraints);

% define time windows
dateGA20 = Leg1.bestGAdate - days(800);
dateGA2End = Leg1.bestGAdate + days(800);

dateArr0 = dateGA20 + 1.5*Leg2.hohmann.TOF_days;
dateArrEnd = dateGA2End + Leg2.hohmann.TOF_days + 5000;

dateSpaceGA2 = linspace(dateGA20,dateGA2End,dateNum);
dateSpaceArr = linspace(dateArr0,dateArrEnd,dateNum);

[GA.V_plus,GA.vInf_plus,...
    Leg2.BestGAdate,Leg2.bestArrdate,GA.matDvArr] = secondLegsynodicAnalysis(constants, data, ...
                                                        dateSpaceGA2,dateSpaceArr, ...
                                                        Leg2.vInf_ga_max, Leg2.vInf_arr_max);

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
plot(date_vec, phase_Leg1, 'b-', 'LineWidth', 1.5)
hold on 
grid on
title('Geometry Repetition: Leg 1 (Jupiter to Mars)')
axis tight
ylabel('Jupiter - Mars [deg]')
ylim([0 360])
datetick('x', 'yyyy')

subplot(2,1,2)
plot(date_vec, phase_Leg2, 'r-', 'LineWidth', 1.5)
hold on
grid on
title('Geometry Repetition: Leg 2 (Mars to Asteroid)')
axis tight
ylabel('Asteroid - Mars Phase [deg]')
ylim([0 360])
datetick('x', 'yyyy')