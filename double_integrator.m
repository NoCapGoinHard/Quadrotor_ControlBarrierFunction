%% November 2025
%% Carlo Rugiero, Francesco Maria Germano, Giovanni Pio Cuoco, Matteo Trusiani

%% Control Barrier Function for collision avoidance motion
%% for an unconstrained robot modeled as a double integrator 

clc; close all; clear variables;

% parameters
global R w mu delta delta1 kp kd next_print_t

next_print_t = 0.01; % printing rate of the debugging string

initialConditions=[1.25;0.25;0;0;0;0]; % [x y z vx vy vz]'

% motion parameter (planar circumpherence centered at the origin)
R = 1; % radious of the circumpherence 
w = 1; % angular velocity

T = 2*pi/w; % simulation length

% reference controller parameters
kp = 100; % proportional gain
kd = 20; % derivative gain

% control barrier function (cbf) parameters
delta = 0.05; % collision thereshold
delta1 = 0.25; % cbf activation thereshold
mu = 0.05; % cbf gain

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

%%

function state_dot = motion_model(t, state)

global mu delta delta1 kp kd next_print_t
p = state(1:3); % robot position
p_dot = state(4:6); % robot velocity

% obstacle
[obs, obs_dot, obs_2dot]=obs_traj(t);

% checks for collisions 
z = p-obs; z_dot = p_dot-obs_dot;
if norm(z)<=delta
    error(['Collision has happened at t = ', num2str(t)]);
end

% trajectory planning
[pd, pd_dot, pd_2dot]=traj_plan(t);

% trajectory tracking controller (reference controller in absence of obstacles)
u_star = pd_2dot + kd*(pd_dot-p_dot) + kp*(pd-p);

% control barrier function (controller with obstacles)
h = z'*(z+mu*z_dot);
proj = (z*z')/(z'*z);
u_cbf = obs_2dot -(2/mu)*proj*z_dot + (eye(3)-proj)*u_star;

% switching condition
% The "if" statement is avoided beacuse it maked the controller discontnuous and ode45 fails.
% The coefficients "a" and "b" generate a numerical AND operator,
% making the control law smooth enough

h_dot_star = z'*(2*z_dot + mu*u_star);
a=min(exp(-100*(h-delta1)),1); % "a" is equivalent to check "h<=delta1": a=1 if h<=delta1, a=0 if h>delta1
b=min(exp(-100*h_dot_star),1); % "b" is equivalent to check "h_dot_star<=0"

% if coeff=0 then u=u_star is active, if coeff=1 then u=u_cbf 
coeff=a*b; 
u = coeff*u_cbf + (1-coeff)*u_star;

% print the current relevant datas for debugging 
if t >= next_print_t
    disp(['t = ', num2str(t), ...
        ', h = ', num2str(h), ...
        ', h_dot_star = ', num2str(h_dot_star), ...
        ', coeff = ', num2str(coeff), ...
        '.']);
    next_print_t = next_print_t + 0.01;
end

% output of the motion model function 
state_dot=[p_dot;u];
end

function [pd, pd_dot, pd_2dot]=traj_plan(t)
% generates a reference trajectory to follow
global R w
pd      = [       R*cos(w*t);        R*sin(w*t); 0];
pd_dot  = [    -R*w*sin(w*t);      R*w*cos(w*t); 0];
pd_2dot = [-R*(w^2)*cos(w*t); -R*(w^2)*sin(w*t); 0];
end

function [obs, obs_dot, obs_2dot]=obs_traj(t)
% generates the obstacle trajectory
obs      = [-1.5*pi+t; -1.15; 0];
obs_dot  = [1; 0; 0];
obs_2dot = [0; 0; 0];
end