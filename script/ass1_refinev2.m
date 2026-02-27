%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- FMINCON OPTIMIZATION--- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial Guess
best_dep = date2mjd2000(datetime2vec(datetime(best_DepDate))); % Data Partenza Giove
best_ga = date2mjd2000(datetime2vec(datetime(best_GADate))); % Data Gravity Assist Marte
best_arr = date2mjd2000(datetime2vec(datetime(best_ArrDate))); % Data Arrivo Asteroide

x0 = [best_dep, best_ga, best_arr];

% Bounds
buffer = 1000; 
lb = x0 - buffer;
ub = x0 + buffer;

% Solver
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...       % Sequential Quadratic Programming è ottimo per traiettorie   
    'StepTolerance', 1e-14, ...
    'ConstraintTolerance', 1e-14, ...
    'MaxFunctionEvaluations', 2000);
%   'Display', 'iter', ...        

% Optimization
fprintf('\nStart Optimization...\n')
tic;

functionCost = @(x) missionCost(x,constants,data);
functionCon = @(x) missionCon(x,constants,data);
[x_opt, fval, exitflag, output] = fmincon(functionCost, x0, [], [], [], [], lb, ub, functionCon, options);
toc;

tDep_opt = datetime(mjd20002date(x_opt(1)));
tGA_opt = datetime(mjd20002date(x_opt(2)));
tArr_opt = datetime(mjd20002date(x_opt(3)));

fprintf('\n--- OPTIMIZED SOLUTION ---\n')
fprintf('Jupiter Departure:         %s\n', tDep_opt)
fprintf('Mars FlyBy:                %s\n', tGA_opt)
fprintf('Asteroid Arrival:          %s\n', tArr_opt)
fprintf('Total Delta v:             %.4f km/s\n', fval)


[~, details] = missionCost(x_opt,constants,data);
fprintf('Departure Delta v:         %.4f km/s\n',details.Dv_Dep)
fprintf('Gravity Assist Delta v:    %.4f km/s\n',details.Dv_GA)
fprintf('Arrival  Delta v:          %.4f km/s\n',details.Dv_Arr)
fprintf('Radius of Pericenter:      %.1f km (Altitude: %.1f km)\n', details.rp, details.rp - 3396);

% results in MJD2000
DEP_mjd2000 = date2mjd2000(datevec(tDep_opt));
GA_mjd2000 = date2mjd2000(datevec(tGA_opt));
ARR_mjd2000 = date2mjd2000(datevec(tArr_opt));
fprintf('--- MJD2000 DATE CONVERSION---\n')
fprintf('Jupiter Departure:         %.4f\n',DEP_mjd2000);
fprintf('Mars Gravity Assist:         %.4f\n',GA_mjd2000);
fprintf('Asteroid Arrival:         %.4f\n',ARR_mjd2000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot3Opt(constants,data,tDep_opt,tGA_opt,tArr_opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GA %%
vInf_minus = details.vInf_in;
vInf_plus = details.vInf_out;
GA = computeGA(data,vInf_minus,vInf_plus);
plotGAOpt(data,GA)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J, info] = missionCost(x,constants,data)
% missionCost Objective function for the  optimization.
%
% DESCRIPTION:
%   Computes the total Delta V cost (Objective Function J) for a mission 
%   composed of two legs with an intermediate Gravity Assist (GA).
%   The function:
%   1. Solves the Lambert problem for Leg 1 (Departure -> Flyby) and Leg 2 
%      (Flyby -> Arrival) based on the epochs in 'x'.
%   2. Evaluates the Gravity Assist geometry at the intermediate planet.
%      - If the required turning angle is feasible naturally, it calculates 
%        the cost as the difference in V_inf magnitude (Ballistic GA).
%      - If the required angle exceeds the physical limit, it calculates the 
%        Delta V required at pericenter to force the turn (Powered GA).
%   3. Sums the impulses: J = DV_Launch + DV_GA + DV_Arrival.
%
% INPUT:
%   x[3x1]          Design variable vector containing epochs [MJD2000]:
%                   [t_Departure; t_Flyby; t_Arrival]
%   constants       Structure with physical constants (muSun, etc.)
%   data            Structure with mission data (IDs, radii, mu, limits)
%
% OUTPUT:
%   J[1]            Total scalar cost (Delta V) [km/s]
%   info[struct]    (Optional) Detailed breakdown of the solution:
%                   .Dv_Dep, .Dv_GA, .Dv_Arr, .rp, .ga_type, etc.
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------
    tDep = x(1); % MJD2000 Dep
    tGA = x(2); % MJD2000 Flyby
    tArr = x(3); % MJD2000 Arr


    [kepDepVec] = uplanet(tDep, data.idDep);
    kepDep = kepVec2Struct(kepDepVec);
    sDep = kep2car(constants.muSun,kepDep);
    rDep = sDep(1:3);
    vDep = sDep(4:6);

    [kepGAVec] = uplanet(tGA, data.idGA);
    kepGA = kepVec2Struct(kepGAVec);
    sGA = kep2car(constants.muSun,kepGA);
    rGA = sGA(1:3);
    vGA = sGA(4:6);

    tofLeg1 = (tGA - tDep) * 86400;

    % Lambert Leg 1
    [~,~,~,~, v1_dep, v1_ga, ~, ~] = lambertMR(rDep, rGA, tofLeg1, constants.muSun, 0, 0, 0);
    
    Dv_Dep = norm(v1_dep' - vDep);
    
    vInf_in = v1_ga' - vGA;
    
    % LEG 2: Mars -> Asteroid 
    kepAstVec = ephAsteroids(tArr, data.asteroidNumber);
    kepAst = kepVec2Struct(kepAstVec);
    sAst = kep2car(constants.muSun,kepAst);
    rAst = sAst(1:3);
    vAst = sAst(4:6);
    
    tofLeg2 = (tArr - tGA) * 86400;

    % Lambert Leg 2
    [~,~,~,~, v2_dep, v2_arr, ~, ~] = lambertMR(rGA, rAst, tofLeg2, constants.muSun, 0, 0, 0);
    
    vInf_out = v2_dep' - vGA;
    
    Dv_Arr = norm(vAst - v2_arr');
    
    % Gravity Assist
    v_in_mag  = norm(vInf_in);
    v_out_mag = norm(vInf_out);

    cos_delta = dot(vInf_in, vInf_out) / (v_in_mag * v_out_mag);
    delta_req = acos(max(-1, min(1, cos_delta)));
    
    r_safe = data.r_Mars + data.h_atmMars; 
    mu = data.muMars;
    
    e_in  = 1 + (r_safe * v_in_mag^2) / mu;
    e_out = 1 + (r_safe * v_out_mag^2) / mu;
    delta_grav = asin(1/e_in) + asin(1/e_out);
    
    if delta_req <= delta_grav
        Dv_GA = abs(v_out_mag - v_in_mag);
        v_mean = (v_in_mag + v_out_mag) / 2;
        e_req = 1 / sin(delta_req/2);
        rp_used = (mu / v_mean^2) * (e_req - 1);
        ga_type = 'Ballistic';
        
    else
        delta_maneuver = delta_req - delta_grav;
        vp_in  = sqrt(v_in_mag^2 + 2*mu/r_safe);
        vp_out = sqrt(v_out_mag^2 + 2*mu/r_safe);        
        Dv_GA = sqrt(vp_in^2 + vp_out^2 - 2*vp_in*vp_out*cos(delta_maneuver));        
        rp_used = r_safe + 50;
        ga_type = 'Powered';
    end

    J = Dv_Dep + Dv_GA + Dv_Arr;
    
    % extra info
    if nargout > 1
        info.Dv_Dep = Dv_Dep;
        info.Dv_GA  = Dv_GA;
        info.Dv_Arr = Dv_Arr;
        info.vInf_in = vInf_in;
        info.vInf_out = vInf_out;
        info.rp = rp_used;
        info.h_flyby = rp_used - data.r_Mars;
        info.ga_type = ga_type;
    end

end


function [c, ceq] = missionCon(x,constants,data)
% missionCon Non-linear constraint function for the optimizer.
%
% DESCRIPTION:
%   Evaluates the non-linear inequality and equality constraints for the 
%   trajectory optimization problem.
%   It calls 'missionCost' to retrieve the computed Delta V values and 
%   enforces operational limits such as:
%   - Maximum characteristic energy (C3) or Delta V at departure.
%   - Maximum Delta V allowed for the arrival braking maneuver.
%
% INPUT:
%   x[3x1]          Design variable vector (Epochs in MJD2000)
%   constants, data Mission configuration structures
%
% OUTPUT:
%   c[Nx1]          Inequality constraints (c <= 0)
%   ceq[Mx1]        Equality constraints (ceq = 0)
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------
    [~, info] = missionCost(x, constants, data);
    
    % Launcher Limit (C3 Max)
    max_v_dep = 15.0; 
    c(1) = info.Dv_Dep - max_v_dep;
    
    % Arrival Limit (Residual Propellant)  
    max_v_arr = 15.0;
    c(2) = info.Dv_Arr - max_v_arr;

    ceq = [];
end

