%% November 2025
%% Carlo Rugiero, Francesco Maria Germano, Giovanni Pio Cuoco, Matteo Trusiani

%% Control Barrier Function for collision avoidance motion
%% for an unconstrained robot modeled as a double integrator

clc; close all; clear variables;

% parameters
global T_s V P_obs obs_est obs_dot_est obs_2dot_est R w mu delta delta1 kp kd traj_type obs_type next_print_t_1 next_print_t_2

% printing rate of the debugging string
next_print_t_1 = 0.1;
next_print_t_2 = 0.1;

% robot initial conditions. state=[x y z vx vy vz]'
initialConditions=[1.25;0.25;0;0;0;0];

traj_type = 'circle';     % 'line' or 'circle'
obs_type = 'circle';      % 'point' or 'line' or 'circle' or 'parabola'

% motion parameter for ref_case=2 (planar circumpherence centered at the origin)
R = 1; % radious of the circumpherence
w = 1; % angular velocity

if strcmp(traj_type,'circle')
    T = 2*pi/w;
else
    T = 15;
end
T_s=0.005; % sampling time

% reference controller parameters
kp = 100; % proportional gain
kd = 20; % derivative gain

% control barrier function (cbf) parameters
delta1 = 0.15; % cbf activation thereshold
delta = delta1/10; % collision thereshold
mu = 0.5; % cbf gain

% Kalman filter parameters
V=0.001*eye(9); % process noise covariance
P_obs=V; % initial value of estimate coviariance
obs_est=obs_traj(0); % initial value of obstacle position estimate
obs_dot_est=zeros(3,1); % initial value of obstacle velocity estimate
obs_2dot_est=zeros(3,1); % initial value of obstacle acceleration estimate

%% running ode
state_current = initialConditions;
t_current = 0;
steps = 1:T/T_s;
t=[]; state=[];

for i=steps
u=controller(t_current,state_current);
[t_ode,state_ode]=ode45(@(t_ode,state_ode) [state_ode(4:6);u], [(i-1)*T_s (i)*T_s], state_current, odeset('RelTol',1e-9,'AbsTol',1e-15));

% update time and state vector
t = [t;t_ode];
state = [state;state_ode];

% update for next iteration 
state_current=state(end,:)';
t_current=t(end);
end

%% get the results
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
    obs=obs_traj(t(i));
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

%% functions

function u=controller(t,state)
global mu delta delta1 kp kd next_print_t_2
p = state(1:3); % robot position
p_dot = state(4:6); % robot velocity

% obstacle
obs=obs_traj(t);
[obs_dot, obs_2dot]=kalman(t,obs);

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
h_dot_star = z'*(2*z_dot + mu*u_star);

if (h<=delta1)&&(h_dot_star<=0)
    u = obs_2dot -(2/mu)*proj*z_dot + (eye(3)-proj)*u_star;
    coeff=1;
else
    u=u_star;
    coeff=0;
end

% print the current relevant datas for debugging
if t >= next_print_t_2
    disp(['t = ', num2str(t), ...
        ', h = ', num2str(h), ...
        ', h_dot_star = ', num2str(h_dot_star), ...
        ', coeff = ', num2str(coeff), ...
        '.']);
    next_print_t_2 = next_print_t_2 + 0.1;
end
end

function [pd, pd_dot, pd_2dot]=traj_plan(t)
% generates a reference trajectory to follow
global traj_type R w
switch traj_type
    case 'line'
        pd      = [t; 0; 0];
        pd_dot  = [1; 0; 0];
        pd_2dot = [0; 0; 0];
    case 'circle'
        pd      = [       R*cos(w*t);        R*sin(w*t); 0];
        pd_dot  = [    -R*w*sin(w*t);      R*w*cos(w*t); 0];
        pd_2dot = [-R*(w^2)*cos(w*t); -R*(w^2)*sin(w*t); 0];
    otherwise
        error('Please select ref_type among the available values');
end
end

function obs=obs_traj(t)
% generates the obstacle trajectory
global obs_type w R
switch obs_type
    case 'point'
        obs      = [7; 0; 0];
    case 'line'
        obs      = [10; 10-t; 0];
    case 'circle'
        r=0.15;
        obs      = [R*(-0.9-r*cos(w*t)); r*R*sin(w*t); 0];
    case 'parabola'
        obs      = [t; t^2-8*t+15; 0];
    otherwise
        error('Please select obs_type among the available values');
end
end

function [obs_dot, obs_2dot]=kalman(t,obs)
global T_s P_obs V obs_est obs_dot_est obs_2dot_est next_print_t_1

% discrete model of the obstacle motion
A = [eye(3), T_s*eye(3), ((T_s^2)/2)*eye(3);
         zeros(3), eye(3), T_s*eye(3);
         zeros(3), zeros(3), eye(3)];
C = [eye(3) zeros([3 6])];

W=0.01*eye(3);  % measurement noise

% prediction step
x_pred=A*[obs_est;obs_dot_est;obs_2dot_est];
P_pred=A*P_obs*A'+V;
Gain=P_pred*C'/(C*P_pred*C'+W); % Kalman gain

% correction step
Inn=obs-C*x_pred; % innovation
x_corr=x_pred+Gain*Inn;
P_obs=(eye(9)-Gain*C)*P_pred;

% update variables for the next iteration
obs_est=x_corr(1:3);
obs_dot_est=x_corr(4:6);
obs_2dot_est=x_corr(7:9);

% print prediction error for debugging
if t >= next_print_t_1
    disp(['|est_error| = ',num2str(norm(obs-obs_est))]);
    next_print_t_1 = next_print_t_1 + 0.1;
end

% output of the function
obs_dot=obs_dot_est;
obs_2dot=obs_2dot_est;
end