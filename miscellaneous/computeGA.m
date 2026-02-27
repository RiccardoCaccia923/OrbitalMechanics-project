function [GA,fval] = computeGA(data,vInf_minus,vInf_plus)
%
% DESCRIPTION:
%   Calculates the geometry and propellant cost of a gravity assist maneuver. 
%
% INPUT:
%   data[struct]        Structure containing planetary and physical constants:
%                       .muMars     [km^3/s^2] Mars gravitational parameter
%                       .r_Mars     [km]       Mars mean radius
%                       .h_atmMars  [km]       Atmosphere height
%                       .muSun      [km^3/s^2] Sun gravitational parameter
%                       .R_Mars     [km]       Mean distance Sun-Mars
%   vInf_minus[3x1]     Incoming hyperbolic excess velocity vector [km/s]
%   vInf_plus[3x1]      Outgoing hyperbolic excess velocity vector [km/s]
%
% OUTPUT:
%   GA[struct]          Structure containing maneuver results:
%                       .rp             [km]   Radius of pericenter
%                       .deltaTot       [rad]  Total required deflection angle
%                       .delta_grav     [rad]  Natural gravity deflection
%                       .delta_maneuver [rad]  Deflection provided by engine
%                       .dv_ratio       [-]    Efficiency ratio (DeltaV/DeltaV_tot)
%                       .T_soi_hours    [h]    Time spent inside the SOI
%                       
%   fval[1]             Residual value from the fzero root finding process [-]
%
% CONTRIBUTORS:
%   Riccardo Caccia
%
% -------------------------------------------------------------------------

    deltaTot = acos(dot(vInf_minus,vInf_plus)/(norm(vInf_minus)*norm(vInf_plus)));
    
    % compute rp with fzero
    rpTarget = @(rp) rpTarget_function(rp,vInf_minus,vInf_plus,data.muMars,deltaTot);
    rp0_guess = data.r_Mars + data.h_atmMars; 
    options = optimset('Display','off');
    [rp_sol, fval, exitflag] = fzero(rpTarget, rp0_guess, options);   
    if exitflag ~= 1
        error('Solver has not found any solution for rp!');
    end
    
    % check on rp
    rp_safe = data.r_Mars + data.h_atmMars;
    rp_limit = rp_safe + 50; 

    if rp_sol <= rp_safe
        % not safe rp -> POWERED GA
        rp_final = rp_limit;
        e_in  = 1 + (rp_final * norm(vInf_minus)^2) / data.muMars;
        e_out = 1 + (rp_final * norm(vInf_plus)^2) / data.muMars;
        delta_grav = asin(1/e_in) + asin(1/e_out);
        delta_maneuver = deltaTot - delta_grav;
        type = 'POWERED';
        
    else
        % safe rp -> BALLISTIC GA
        rp_final = rp_sol;
        delta_maneuver = 0; 
        type = 'BALLISTIC';
    end

    vp_minus = sqrt(norm(vInf_minus)^2 + 2*data.muMars/rp_final); 
    vp_plus  = sqrt(norm(vInf_plus)^2 + 2*data.muMars/rp_final); 
    
    delta_vp = sqrt(vp_minus^2 + vp_plus^2 - 2*vp_minus*vp_plus*cos(delta_maneuver));
    
    DvTot_vec = vInf_plus - vInf_minus;
    DvTot = norm(DvTot_vec);
    
    gain = DvTot - delta_vp;
    dv_ratio = delta_vp / DvTot;

    r_soi = data.R_Mars * (data.muMars / data.muSun)^(2/5);
    t_in  = hyperbolicTimeLaw(norm(vInf_minus), rp_final, data.muMars, r_soi);
    t_out = hyperbolicTimeLaw(norm(vInf_plus), rp_final, data.muMars, r_soi);
    T_soi_sec = t_in + t_out;
    T_soi_hours = T_soi_sec / 3600;
    T_soi_days = T_soi_sec / 86400;

    
    fprintf('\n--- GA MANOEUVER RESULTS (%s) ---\n', type)
    fprintf('Radius of Perigee:         %.3f km (Safe: %.1f)\n', rp_final, rp_limit)
    fprintf('Height of passage:         %.3f km\n', rp_final - data.r_Mars);
    fprintf('Time inside SOI:           %.3f hours\n', T_soi_hours);
    fprintf('Total Turn Angle Req:      %.2f deg\n', rad2deg(deltaTot))
    if delta_maneuver > 0
        fprintf('Gravity Turn:              %.2f deg\n', rad2deg(delta_grav))
        fprintf('Engine Turn Needed:        %.2f deg\n', rad2deg(delta_maneuver))
    end
    fprintf('Vp IN:                     %.3f km/s\n', vp_minus) 
    fprintf('Vp OUT:                    %.3f km/s\n', vp_plus) 
    fprintf('Delta V "Paid" (Engine):   %.3f km/s\n', delta_vp)
    fprintf('Delta V "Obtained" (Tot):  %.3f km/s\n', DvTot)
    fprintf('---------------------------------------\n')
    fprintf('Impulse/Total Ratio:       %.4f\n', dv_ratio);
    fprintf('Delta V GAIN:              %.3f km/s\n', gain)

    GA.rp = rp_final;
    GA.deltaTot = deltaTot;
    GA.delta_grav = delta_grav;
    GA.delta_maneuver = delta_maneuver;
    GA.vInf_minus = vInf_minus;
    GA.vInf_plus = vInf_plus;
    GA.vp_minus = vp_minus;
    GA.vp_plus = vp_plus;
    GA.T_soi_hours = T_soi_hours; 
    GA.dv_ratio = dv_ratio;

end