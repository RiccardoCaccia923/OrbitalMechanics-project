function plotOribtParametersV2(t,values,tPert,valuesPert)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION:
%   Generates a comprehensive set of plots to visualize the time evolution 
%   of orbital elements and constants of motion.
%   - Default Mode (nargin != 4): Plots angular momentum (h), eccentricity (e), 
%     specific energy, and velocity components to verify conservation laws 
%     in unperturbed motion.
%   - Comparison Mode (nargin == 4): Overlays perturbed vs. unperturbed 
%     data to analyze the impact of perturbations on the orbital geometry 
%     and to visualize the variation of osculating elements over time.
%
% INPUT
%   t: integration timestamp  
%   s: [r,v] 6x1 state vector
%   values: struct with fields
%           -   type:   type of orbit               'string'
%           -   h:      angular momentum            Nx3 matrix
%           -   hNorm:  norm of angular momentum    Nx1 vector
%           -   e:      eccentricity vector         Nx3 matrix
%           -   eNorm:  eccentricity vector norm    Nx1 vector
%           -   err_he: misOrtogonality of h and e  Nx1 vector
%           -   eps:    specific mechanical energy  Nx1 vector
%           -   v_r:    radial velocity             Nx1 vector
%           -   v_t:    transversal velocity        Nx1 vector
%           -   a:      semi-major axis             scalar
%           -   p:      semi-latus rectum           scalar
%           -   T:      orbital period              scalar
%
% CONTRIBUTORS:
%   Riccardo Caccia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 4
    % plotOrbitParametersV1 --> only unperturbed 
    figure
    subplot(3,1,1)
    plot(t/values.T,values.hNorm)
    grid on
    title('Angular Momentum norm conservation')
    xlabel('t/T')
    ylabel('|h|')
    
    subplot(3,1,2)
    plot(t/values.T,values.eNorm)
    grid on
    title('Eccentricity norm conservation')
    xlabel('t/T')
    ylabel('|e|')
    %ylim([-1,1])
    
    subplot(3,1,3)
    plot(t/values.T,values.err_he)
    grid on
    title('Error for h and e to be ortogonal')
    xlabel('t/T')
    ylabel('error')
    
    
    figure
    subplot(1,2,1)
    plot(t/values.T,values.hNorm,'k--','LineWidth',2)
    hold on
    plot(t/values.T,values.h(1,:),'b--','LineWidth',2)
    plot(t/values.T,values.h(2,:),'r--','LineWidth',2)
    plot(t/values.T,values.h(3,:),'g--','LineWidth',2)
    grid on
    title('Angular Momentum comp. conservation')
    xlabel('t/T', 'Interpreter', 'latex')
    ylabel('$|h|$ , $h_x$ , $h_y$ , $h_z$', 'Interpreter', 'latex')
    legend('$|h|$' , '$h_x$' , '$h_y$' , '$h_z$', 'Interpreter', 'latex')
    
    subplot(1,2,2)
    plot(t/values.T,values.eNorm,'k--','LineWidth',2)
    hold on
    plot(t/values.T,values.e(1,:),'b:','LineWidth',2)
    plot(t/values.T,values.e(2,:),'r--','LineWidth',2)
    plot(t/values.T,values.e(3,:),'g:','LineWidth',2)
    grid on
    title('Eccentricity comp. conservation')
    xlabel('t/T', 'Interpreter', 'latex')
    ylabel('$|e|$ , $e_x$ , $e_y$ , $e_z$', 'Interpreter', 'latex')
    legend('$|e|$' , '$e_x$' , '$e_y$' , '$e_z$', 'Interpreter', 'latex')
    
    
    figure
    plot(t/values.T,values.eps)
    grid on
    title('Specific Mechanical Energy conservation')
    xlabel('t/T')
    ylabel('eps[]')
    
    
    figure
    subplot(3,1,1);
    plot(t/values.T,values.v_r)
    hold on
    plot(t/values.T,values.v_t)
    title('Radial & Transversal Velocity')
    xlabel('t/T')
    ylabel('$v_r [m/s]$ , $v_t [m/s]$', 'Interpreter', 'latex')
    legend('$v_r [m/s]$', '$v_t [m/s]$', 'Interpreter', 'latex')
    
    subplot(3,1,2);
    plot(t/values.T,values.v_r)
    grid on
    title('Radial Velocity')
    xlabel('t/T')
    ylabel('$v_r [m/s]$', 'Interpreter', 'latex')
    
    subplot(3,1,3);
    plot(t/values.T,values.v_t)
    grid on
    title('Transversal Velocity')
    xlabel('t/T')
    ylabel('$v_t [m/s]$', 'Interpreter', 'latex')
    
else
    
    % plot angular momentum and eccentricity norm 
    figure
    subplot(2,1,1)
    plot(t/values.T,values.hNorm)
    hold on
    plot(tPert/valuesPert.T,valuesPert.hNorm)
    grid on
    title('Angular Momentum norm conservation')
    xlabel('t/T')
    ylabel('$|h|$', 'Interpreter', 'latex')
    legend('$|h|$','$|h_{pert}|$', 'Interpreter', 'latex')

    subplot(2,1,2)
    plot(t/values.T,values.eNorm)
    hold on
    plot(tPert/valuesPert.T,valuesPert.eNorm)
    grid on
    title('Eccentricity norm conservation')
    xlabel('t/T')
    ylabel('$|e|$', 'Interpreter', 'latex')
    legend('$|e|$','$|e_{pert}|$', 'Interpreter', 'latex')


    figure
    plot(t/values.T,values.err_he)
    hold on
    plot(tPert/valuesPert.T,valuesPert.err_he)
    grid on
    title('Error for h and e to be ortogonal')
    xlabel('t/T','Interpreter', 'latex')
    ylabel('error','Interpreter', 'latex')
    legend('$err-eh$','$err-eh_{pert}$','Interpreter', 'latex')


    figure
    plot(t/values.T,values.hNorm,'k--','LineWidth',2)
    hold on
    plot(t/values.T,values.h(1,:),'b--','LineWidth',2)
    plot(t/values.T,values.h(2,:),'r--','LineWidth',2)
    plot(t/values.T,values.h(3,:),'g--','LineWidth',2)
    plot(tPert/valuesPert.T,valuesPert.hNorm,'k--','LineWidth',2)
    plot(tPert/valuesPert.T,valuesPert.h(1,:),'b--','LineWidth',2)
    plot(tPert/valuesPert.T,valuesPert.h(2,:),'r--','LineWidth',2)
    plot(tPert/valuesPert.T,valuesPert.h(3,:),'g--','LineWidth',2)
    grid on
    title('Angular Momentum comp. conservation')
    xlabel('t/T', 'Interpreter', 'latex')
    %ylabel(['$|h|$ , $h_x$ , $h_y$ , $h_z$ , $|h_{pert}|$ , $h_{x,pert}$ , ' ...
    %    '$h_{y,pert}$ , $h_{z,pert}$'] , 'Interpreter', 'latex')
    ylabel('$h$ , $h_{pert}$', 'Interpreter', 'latex')
    legend('$|h|$' , '$h_x$' , '$h_y$' , '$h_z$', '$|h_{pert}|$' , ...
        '$h_{x,pert}$' , '$h_{y,pert}$' , '$h_{z,pert}$' , 'Interpreter', 'latex')


    figure
    plot(t/values.T,values.eNorm,'k--','LineWidth',2)
    hold on
    plot(t/values.T,values.e(1,:),'b:','LineWidth',2)
    plot(t/values.T,values.e(2,:),'r--','LineWidth',2)
    plot(t/values.T,values.e(3,:),'g:','LineWidth',2)
    plot(tPert/valuesPert.T,valuesPert.eNorm,'k--','LineWidth',2)
    plot(tPert/valuesPert.T,valuesPert.e(1,:),'b--','LineWidth',2)
    plot(tPert/valuesPert.T,valuesPert.e(2,:),'r--','LineWidth',2)
    plot(tPert/valuesPert.T,valuesPert.e(3,:),'g--','LineWidth',2)
    grid on
    title('Eccentricity comp. conservation')
    xlabel('t/T', 'Interpreter', 'latex')
    % ylabel(['$|e|$ , $e_x$ , $e_y$ , $e_z$ , $|e_{pert}|$ , $e_{x,pert}$ , ' ...
    %     '$e_{y,pert}$ , $e_{z,pert}$'] , 'Interpreter', 'latex')
    ylabel('$e$ , $e_{pert}$', 'Interpreter', 'latex')
    legend('$|e|$' , '$e_x$' , '$e_y$' , '$e_z$', '$|e_{pert}|$' , ...
        '$e_{x,pert}$' , '$e_{y,pert}$' , '$e_{z,pert}$' , 'Interpreter', 'latex')

    figure
    plot(t/values.T,values.eps)
    hold on
    plot(tPert/valuesPert.T,valuesPert.eps)
    grid on
    title('Specific Mechanical Energy conservation')
    xlabel('t/T')
    ylabel('$eps$','Interpreter', 'latex')
    legend('$eps$','$eps_{pert}$','Interpreter', 'latex')

    
    figure
    plot(t/values.T,values.v_r,'g')
    hold on
    plot(t/values.T,values.v_t,'k')
    plot(tPert/valuesPert.T,valuesPert.v_r,'b--')
    plot(tPert/valuesPert.T,valuesPert.v_t,'r--')
    title('$Radial & Transversal Velocity$')
    xlabel('t/T')
    ylabel('$v_r $ , $v_t$ , $v_{r,pert} $ , $v_{t,pert}$ $[m/s]$', ...
        'Interpreter', 'latex')
    legend('$v_r $', '$v_t $', '$v_{r,pert}$' , '$v_{t,pert}$', ...
        'Interpreter', 'latex')

    
    figure
    subplot(2,1,1);
    plot(t/values.T,values.v_r,'g')
    hold on
    plot(tPert/valuesPert.T,valuesPert.v_r,'b--')
    grid on
    title('Radial Velocity')
    xlabel('t/T')
    ylabel('$v_r$ , $v_{r,pert}$ $[m/s]$', 'Interpreter', 'latex')
    legend('$v_r$' , '$v_{r,pert}$', 'Interpreter', 'latex')

    subplot(2,1,2);
    plot(t/values.T,values.v_t,'k')
    hold on
    plot(tPert/valuesPert.T,valuesPert.v_t,'r--')
    grid on
    title('Transversal Velocity')
    xlabel('t/T')
    ylabel('$v_t$ , $v_{t,pert}$ $[m/s]$', 'Interpreter', 'latex')
    legend('$v_t$' , '$v_{t,pert}$', 'Interpreter', 'latex')
end

end 