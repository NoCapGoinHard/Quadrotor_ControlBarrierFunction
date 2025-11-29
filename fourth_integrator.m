%% November 2025
%% Carlo Rugiero, Francesco Maria Germano, Giovanni Pio Cuoco, Matteo Trusiani

%% Control Barrier Function for collision avoidance motion
%% for an unconstrained robot modeled as a double integrator

clc; close all; clear variables;

% parameters
global R w obs_case mu2 mu3 mu4 delta delta1 k0 k1 k2 k3 next_print_t

next_print_t = 0.01; % printing rate of the debugging string

% robot initial conditions. state=[x y z vx vy vz ax ay az jx jy jz]'
initialConditions=[0;0;0; 0;0;0; 0;0;0; 0;0;0];

% obstacle trajectory.
% obs_case=1: fixed on the point (-1,0); obs_case=2: moves on the vertical
% line x=-1; obs_case=3: moves on a cirumpherence near (-1,0)
obs_case = 1;

% motion parameter (planar circumpherence centered at the origin)
R = 1; % radious of the circumpherence
w = 1; % angular velocity

T = 9.2;%2*pi/w; % simulation length

% reference controller parameters
eig=10;
k0 = 1*eig^4;
k1 = 4*eig^3;
k2 = 6*eig^2;
k3 = 4*eig;

% control barrier function (cbf) parameters
delta1 = 1.5; % cbf activation threshold
delta = delta1/10; % collision threshold
mu=0.35;
q=8;
mu2 = mu*q; 
mu3 = 0;
mu4 = mu;

% running ode
[t,state]=ode45(@(t,state) motion_model(t, state), [0 T], initialConditions, odeset('RelTol',1e-9,'AbsTol',1e-15));

% get the results
x = state(:,1);
y = state(:,2);
z = state(:,3);
xd = zeros(length(t),1);
yd = zeros(length(t),1);
%zd = zeros(length(t),1);
x_obs = zeros(length(t),1);
y_obs = zeros(length(t),1);
%y_obs = zeros(length(t),1);

for i = 1:length(t)
    [pd, ~, ~]=traj_plan(t(i));
    xd(i)   = pd(1);
    yd(i)   = pd(2);
    %    zd(i)   = pd(3);
    [obs, ~, ~]=obs_traj(t(i));
    x_obs(i) = obs(1);
    y_obs(i) = obs(2);
    %    z_obs = obs(3);
end

% plot planar trajectory on the x-y plane
figure(1);
image = plot(x, y, xd, yd, x_obs, y_obs);
set(image(1), 'LineStyle', '-', 'Color', [0 0 0.55], 'LineWidth', 1.5);
set(image(2), 'LineStyle', '--','Color', [1 0.4 0.2], 'LineWidth', 1.5);
set(image(3), 'LineStyle', '--','Color', [0 0 0], 'LineWidth', 5);
legend('y','y_d','obstacle');
xlabel('x [m]'), ylabel('y [m]');
title('Position: y(x) and yd(xd)');
figure(2);
plot(t, x);

%%

function state_dot = motion_model(t, state)

global mu2 mu3 mu4 delta delta1 k0 k1 k2 k3 next_print_t
p = state(1:3); % robot position
p_dot = state(4:6); % robot velocity
p_2dot = state(7:9); % robot acceleration
p_3dot = state(10:12); % robot jerk

% obstacle
[obs, obs_dot, obs_2dot, obs_3dot, obs_4dot]=obs_traj(t);

% checks for collisions
z = p-obs; z_dot = p_dot-obs_dot; z_2dot = p_2dot-obs_2dot; z_3dot = p_3dot-obs_3dot;
if norm(z)<=delta
    error(['Collision has happened at t = ', num2str(t)]);
end

% trajectory planning
[pd, pd_dot, pd_2dot, pd_3dot, pd_4dot]=traj_plan(t);

% trajectory tracking controller (reference controller in absence of obstacles)
u_star = pd_4dot + k3*(pd_3dot-p_3dot) + k2*(pd_2dot-p_2dot) + k1*(pd_dot-p_dot) + k0*(pd-p);

% control barrier function (controller with obstacles)
h=z'*(z+mu2*z_dot+mu3*z_2dot+mu4*z_3dot);
u_cbf = obs_4dot -(2*z_dot+mu2*z_2dot+mu3*z_3dot)/mu4 -(mu3*z_dot'*z_2dot+mu4*z_dot'*z_3dot)*(z+mu2*z_dot+mu3*z_2dot+mu4*z_3dot)/(mu4*h) + (eye(3)-(z*((z.'*z)^-1)*z.'))*u_star;

% switching condition
% The "if" statement is avoided beacuse it maked the controller discontnuous and ode45 fails.
% The coefficients "a" and "b" generate a numerical AND operator,
% making the control law smooth enough

h_dot_star=z_dot'*(z+mu2*z_dot+mu3*z_2dot+mu4*z_3dot)+z'*(z_dot+mu2*z_2dot+mu3*z_3dot+mu4*(u_star-obs_4dot));
a=min(exp(-100*(h-delta1)),1); % "a" is equivalent to check "h<=delta1": a=1 if h<=delta1, a=0 if h>delta1
b=min(exp(-100*h_dot_star),1); % "b" is equivalent to check "h_dot_star<=0"

% if coeff=0 then u=u_star is active, if coeff=1 then u=u_cbf
coeff=a*b;
u = coeff*u_cbf + (1-coeff)*u_star;

% print the current relevant datas for debugging
if t >= next_print_t
    disp(['t = ', num2str(t), ...
        ', norm(z) = ', num2str(norm(z)), ...
        ', h = ', num2str(h), ...
        ', h_dot_star = ', num2str(h_dot_star), ...
        ', coeff = ', num2str(coeff), ...
        '.']);
    next_print_t = next_print_t + 0.01;
end

% output of the motion model function
state_dot=[p_dot;p_2dot;p_3dot;u];
end

function [pd, pd_dot, pd_2dot, pd_3dot, pd_4dot]=traj_plan(t)
% generates a reference trajectory to follow
global R w
% pd      = [       R*cos(w*t);        R*sin(w*t); 0];
% pd_dot  = [    -R*w*sin(w*t);      R*w*cos(w*t); 0];
% pd_2dot = [-R*(w^2)*cos(w*t); -R*(w^2)*sin(w*t); 0];
pd      = [t; 0; 0];
pd_dot  = [1; 0; 0];
pd_2dot = [0; 0; 0];
pd_3dot = [0; 0; 0];
pd_4dot = [0; 0; 0];
end

function [obs, obs_dot, obs_2dot, obs_3dot, obs_4dot]=obs_traj(t)
% generates the obstacle trajectory
global obs_case w R

switch obs_case
    case 1 % point
        obs      = [3; 0; 0];
        obs_dot  = [0; 0; 0];
        obs_2dot = [0; 0; 0];
        obs_3dot = [0; 0; 0];
        obs_4dot = [0; 0; 0];
    case 2 % line
        obs      = [-1; pi-t; 0];
        obs_dot  = [0; -1; 0];
        obs_2dot = [0; 0; 0];
        obs_3dot = [0; 0; 0];
        obs_4dot = [0; 0; 0];
    case 3 % circumpherence
        obs      = [R*(-0.9-0.3*cos(w*t)); 0.3*R*sin(w*t); 0];
        obs_dot  = [0.3*R*w*sin(w*t);      0.3*R*w*cos(w*t); 0];
        obs_2dot = [0.3*R*(w^2)*cos(w*t); -0.3*R*(w^2)*sin(w*t); 0];

    otherwise
        error('Please select obs_case among the available values');
end

end