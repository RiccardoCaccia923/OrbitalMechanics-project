function [t,s,values] = twoBodyProblemNoPlot(mu,s0,nPeriods,time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% twoBodyProblem: function created to solve the two bod problem for the
%                 non-perturbed case by numerical integration
% DESCRIPTION:
%   Solves the unperturbed Two-Body Problem via numerical integration.
%   The function propagates the initial state vector [r; v] forward in time
%   using the ode solver (ode113).
%   SPECIFIC VERSION OF twoBodyProblem WITHOUT PLOT ONLY TO AVOID PLOT
%   REDUNDANCY IN SOME ANALYSIS 
%
% INPUTS
%   mu:         gravitational parameter [km^3/s^2]
%   s0:         [r0; v0] initial state vector
%   nPeriods:   number of orbital periods on which perform the analysis
%   time:       timespan not in form of nPeriods
%
% OUTPUT
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

%% extract fields r0, v0
[x,y]=size(s0);
if x ~= 6 && y ~= 1
    disp('ERROR: input s0 not in the right form, need a column vector [6x1]')
end

s0 = s0(:);                 % force s0 in column
r0 = s0(1:3);               % 3x1
v0 = s0(4:6);               % 3x1

r0Norm = norm(r0);
v0Norm = norm(v0);

%% checks
% if nargin ~= 4
%     nPeriodsPar = nPeriodsPlot;
% end

% check r0 concistency
if norm(r0) == 0
    disp('error: r0 is not consistent')
end

% check specific energy to define the type of orbit
eps0 = v0Norm^2/2 - mu/r0Norm;
if eps0 < 0 
    type = 'closed';
elseif eps0 == 0
    type = 'parabolic';
else
    type = 'hyperbolic';
end

%% compute geometric parameters
% compute angular momentum
h0 = cross(r0,v0);                  % vector
h0Norm = norm(h0);                  % scalar

% define eccentricity vector
e0 = (1/mu)*cross(v0,h0)-r0/r0Norm;   % vector
e0Norm = norm(e0);                    % scalar

% compute semilatus rectum 
p = h0Norm^2/mu;                    % scalar

% enter the right case
if strcmp(type, 'closed')

    a = p/(1-e0Norm^2);              % Semi-major axis [km]
    
    T = 2*pi*sqrt(a^3/mu);          % Orbital period [s]

    if isempty(nPeriods)
        nPeriods = time/T;
        % tspan = linspace(0 , time , 1e5);
        tspan = time;
    elseif isempty(time)
        tspan = linspace(0 , nPeriods * T , 1e5);
    end

end

% define state derivatives
twoBody = @(t,s) [s(4:6) ;
                 -mu*s(1:3)/(norm(s(1:3))^3)];

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[t,s] = ode113(twoBody, tspan, s0, options );

%% plot orbit
% x = s(:,1);
% y = s(:,2);
% z = s(:,3);
% 
% 
% plot3(x,y,z,'LineWidth',1.5)
% hold on
% plot3(x(1), y(1), z(1), 'go','MarkerFaceColor','g'); 
% plot3(x(end), y(end), z(end), 'ro','MarkerFaceColor','none','MarkerEdgeColor','r'); 
% 
% plotEarth()

%% compute parameters
r = [s(:,1),s(:,2),s(:,3)]';
v = [s(:,4),s(:,5),s(:,6)]';

h=cross(r,v);
hNorm = zeros(1,length(t));

e = zeros(3,length(t));
eNorm = zeros(1,length(t));

eps = zeros(1,length(t));

v_r = zeros(1,length(t));
v_t = zeros(1,length(t));

 for  i = 1 : length(t)

    hNorm(i) = norm(h(:,i));

    e(:,i) = (1/mu) * cross(v(:,i),h(:,i)) - r(:,i)/norm(r(:,i));
    eNorm(i) = norm(e(:,i));

    eps(i) = 0.5 * norm(v(:,i))^2 - mu/norm(r(:,i));

    u_r = r(:,i) / norm(r(:,i));
    u_h = h(:,i) / hNorm(i);
    u_t = cross(u_h,u_r);

    v_r(i) = dot(v(:,i),u_r);
    v_t(i) = dot(v(:,i),u_t);

 end

% check if h and e are ortogonal
err_he = dot(h,e);


%% create values struct 
values.type = type;
values.h = h;
values.hNorm = hNorm;
values.e = e;
values.eNorm = eNorm;
values.err_he = err_he;
values.eps = eps;
values.v_r = v_r;
values.v_t = v_t;
values.a = a;
values.p = p;
values.T = T;

end