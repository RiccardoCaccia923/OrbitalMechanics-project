function plotOribtParameters(t,values)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           GO TO plotOrbitParametersV2
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
figure
subplot(3,1,1)
plot(t,values.hNorm)
grid on
title('Angular Momentum norm conservation')
xlabel('t[s]')
ylabel('|h|')

subplot(3,1,2)
plot(t,values.eNorm)
grid on
title('Eccentricity norm conservation')
xlabel('t[s]')
ylabel('|e|')
%ylim([-1,1])

subplot(3,1,3)
plot(t,values.err_he)
grid on
title('Error for h and e to be ortogonal')
xlabel('t[s]')
ylabel('error')


figure
subplot(1,2,1)
plot(t,values.hNorm,'k--','LineWidth',2)
hold on
plot(t,values.h(1,:),'b--','LineWidth',2)
plot(t,values.h(2,:),'r--','LineWidth',2)
plot(t,values.h(3,:),'g--','LineWidth',2)
grid on
title('Angular Momentum comp. conservation')
xlabel('t[s]', 'Interpreter', 'latex')
ylabel('$|h|$ , $h_x$ , $h_y$ , $h_z$', 'Interpreter', 'latex')
legend('$|h|$' , '$h_x$' , '$h_y$' , '$h_z$', 'Interpreter', 'latex')

subplot(1,2,2)
plot(t,values.eNorm,'k--','LineWidth',2)
hold on
plot(t,values.e(1,:),'b:','LineWidth',2)
plot(t,values.e(2,:),'r--','LineWidth',2)
plot(t,values.e(3,:),'g:','LineWidth',2)
grid on
title('Eccentricity comp. conservation')
xlabel('t[s]', 'Interpreter', 'latex')
ylabel('$|e|$ , $e_x$ , $e_y$ , $e_z$', 'Interpreter', 'latex')
legend('$|e|$' , '$e_x$' , '$e_y$' , '$e_z$', 'Interpreter', 'latex')


figure
plot(t,values.eps)
grid on
title('Specific Mechanical Energy conservation')
xlabel('t[s]')
ylabel('eps[]')


figure
subplot(3,1,1);
plot(t,values.v_r)
hold on
plot(t,values.v_t)
title('Radial & Transversal Velocity')
xlabel('t[s]')
ylabel('$v_r [m/s]$ , $v_t [m/s]$', 'Interpreter', 'latex')
legend('$v_r [m/s]$', '$v_t [m/s]$', 'Interpreter', 'latex')

subplot(3,1,2);
plot(t,values.v_r)
grid on
title('Radial Velocity')
xlabel('t[s]')
ylabel('$v_r [m/s]$', 'Interpreter', 'latex')

subplot(3,1,3);
plot(t,values.v_t)
grid on
title('Transversal Velocity')
xlabel('t[s]')
ylabel('$v_t [m/s]$', 'Interpreter', 'latex')

end