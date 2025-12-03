%% December 2025
%% Carlo Rugiero, Francesco Maria Germano, Giovanni Pio Cuoco, Matteo Trusiani

%% Control Barrier Function for collision avoidance motion
%% for an unconstrained robot modeled as a double integrator + multiple obstacles

clc; close all; clear variables;

% GLOBAL PARAMETERS
global R w mu delta delta1 kp kd next_print_t M traj_type

traj_type = 'line';     % 'line' or 'circle'

next_print_t = 0.01;     % printing rate of the debugging string

% robot initial conditions. state=[x y z vx vy vz]'
initialConditions = [0; 0; 0; 0; 0; 0];

% tracking trajectory (circle)
R = 1;     
w = 1;     

if strcmp(traj_type,'circle')
    T = 2*pi/w;
else
    T = 10;
end

% reference controller parameters
kp = 100;
kd = 20;

% CBF parameters
delta1 = 0.15;       % CBF activation threshold
delta  = delta1/10;  % collision threshold
mu     = 0.5;

% number of obstacles
M = 3;  

% check condition to avoid overlapping 
[obs0, ~, ~] = obs_traj_multi(0);

for i = 1:M
    for j = i+1:M
        if norm(obs0(:,i) - obs0(:,j))^2 <= 2*delta1
            error('Obstacles violate condition ||p_i - p_j||^2 > 2 delta1');
        end
    end
end


%running ode
[t,state] = ode45(@(t,state) motion_model_multi(t, state),[0 T], initialConditions,odeset('RelTol',1e-9,'AbsTol',1e-15));

% get the results
x = state(:,1);
y = state(:,2);

xd   = zeros(length(t),1);
yd   = zeros(length(t),1);

xobs = zeros(length(t),M);
yobs = zeros(length(t),M);

d_min  = zeros(length(t),1);

% compute reference and obstacles trajectories
for i = 1:length(t)

    % reference
    [pd, ~, ~] = traj_plan(t(i));
    xd(i) = pd(1);
    yd(i) = pd(2);

    % obstacles
    [obs_all, ~, ~] = obs_traj_multi(t(i));

    dist2 = zeros(1,M);

    for j = 1:M
        xobs(i,j) = obs_all(1,j); % the x-coordinate of obstacle j at time t(i)
        yobs(i,j) = obs_all(2,j); % the y-coordinate of obstacle j at time t(i)
        dist2(j)  = (state(i,1) - obs_all(1,j))^2 + (state(i,2) - obs_all(2,j))^2; % squared Euclidean distance = (x_robot - x_obstacle)^2 + (y_robot - y_obstacle)^2
    end

    % minimum distance
    d_min(i) = sqrt(min(dist2));
end


figure(1); hold on; grid on; axis equal;
plot(x, y, 'Color', [0 0 0.55], 'LineWidth', 1.8); 
plot(xd, yd, '--', 'Color', [1 0.4 0.2], 'LineWidth', 1.5);

for j = 1:M
    plot(xobs(:,j), yobs(:,j), '--','Color',[0 0 0], 'LineWidth', 2);
end

legend('Robot', 'Reference', 'Obstacles');
xlabel('x [m]'); ylabel('y [m]');
title(['Robot trajectory with multiple obstacles']);

% plot minum distance vs time
figure(2); hold on; grid on;

plot(t, d_min, 'b', 'LineWidth', 1.8);

yline(delta,  'r--', 'LineWidth', 2);
yline(delta1, 'g--', 'LineWidth', 2);

xlabel('time [s]');
ylabel('min distance [m]');

legend('minimum distance', 'collision threshold \delta', 'CBF threshold \delta_1');

title('Minimum distance to obstacles over time');



function state_dot = motion_model_multi(t, state)

global mu delta delta1 kp kd next_print_t M

p     = state(1:3);   % robot position
p_dot = state(4:6);   % robot velocity

% get all the obstacles
[obs_all, obs_dot_all, obs_2dot_all] = obs_traj_multi(t);

% compute relative vectors and distances
z_all     = zeros(3,M);
z_dot_all = zeros(3,M);
dist2     = zeros(1,M);

for i = 1:M
    z_all(:,i)     = p     - obs_all(:,i);
    z_dot_all(:,i) = p_dot - obs_dot_all(:,i); 
    dist2(i)       = z_all(:,i)'*z_all(:,i); % = ||p - p_obs,i||^2
end

% select the closest obstacle
[~, i_star] = min(dist2);

% active obstacle
obs      = obs_all(:,i_star);
obs_dot  = obs_dot_all(:,i_star);
obs_2dot = obs_2dot_all(:,i_star);

% active relative quantities
z     = z_all(:,i_star);   % z = p-p_obs,i
z_dot = z_dot_all(:,i_star);

% check collision of the active obstacle
if norm(z) <= delta
    error(['Collision happened at t = ', num2str(t)]);
end

[pd, pd_dot, pd_2dot] = traj_plan(t); % trajectory planning

u_star = pd_2dot + kd*(pd_dot-p_dot) + kp*(pd-p);  %trajectory tracking controller (reference controller in absence of obstacles)

% control barrier function (controller with obstacles)
h = z'*(z + mu*z_dot);
proj = (z*z')/(z'*z);

% CBF avoidance control
u_cbf = obs_2dot -(2/mu)*proj*z_dot + (eye(3)-proj)*u_star;

h_dot_star = z'*(2*z_dot + mu*u_star);

a = min(exp(-100*(h-delta1)),1);   % a=1 if h <= delta1
b = min(exp(-100*h_dot_star),1);   % b=1 if h_dot_star <= 0

% if coeff=0 then u=u_star is active, if coeff=1 then u=u_cbf
coeff = a*b;
u = coeff*u_cbf + (1-coeff)*u_star;

% debug 
if t >= next_print_t
    disp(['t = ', num2str(t), ...
        ', h = ', num2str(h), ...
        ', coeff = ', num2str(coeff), ...
        ', active obs = ', num2str(i_star), ...
        ', norm(z) = ', num2str(norm(z)), ...
        '.']);
    next_print_t = next_print_t + 0.01;
end

% output of the motion model function
state_dot = [p_dot; u];

end


% generates a reference trajectory to follow (circle or line)
function [pd, pd_dot, pd_2dot] = traj_plan(t)

global R w traj_type

switch traj_type

    case 'circle'
        pd      = [       R*cos(w*t);        R*sin(w*t); 0];
        pd_dot  = [    -R*w*sin(w*t);      R*w*cos(w*t); 0];
        pd_2dot = [-R*(w^2)*cos(w*t); -R*(w^2)*sin(w*t); 0];

    case 'line'
        pd      = [t; 0; 0];
        pd_dot  = [1; 0; 0];
        pd_2dot = [0; 0; 0];

end

end

% generate multiple obstacle trajectory 
function [obs_all, obs_dot_all, obs_2dot_all] = obs_traj_multi(t)

global w R M delta1 traj_type

obs_all     = zeros(3,M);
obs_dot_all = zeros(3,M);
obs_2dot_all= zeros(3,M);

switch traj_type

    case 'circle'
        % OBSTACLE 1: static
        obs_all(:,1)     = [0; -1; 0];
        obs_dot_all(:,1) = [0; 0; 0];
        obs_2dot_all(:,1)= [0; 0; 0];

        % OBSTACLE 2: vertical line
        obs_all(:,2)     = [-1; pi - t; 0];
        obs_dot_all(:,2) = [0; -1; 0];
        obs_2dot_all(:,2)= [0; 0; 0];

        % OBSTACLE 3: circle 
        r = 0.3;
        obs_all(:,3)     = [R*(0 - r*cos(w*t)); 1 + r*R*sin(w*t); 0];
        obs_dot_all(:,3) = [r*R*w*sin(w*t);       r*R*w*cos(w*t); 0];
        obs_2dot_all(:,3)= [r*R*(w^2)*cos(w*t);  -r*R*(w^2)*sin(w*t); 0];

    case 'line'
        % OBSTACLE 1: static
        obs_all(:,1)     = [4; 0; 0];
        obs_dot_all(:,1) = [0; 0; 0];
        obs_2dot_all(:,1)= [0; 0; 0];

        % OBSTACLE 2: vertical line
        obs_all(:,2)     = [3; pi - t; 0];
        obs_dot_all(:,2) = [0; -1; 0];
        obs_2dot_all(:,2)= [0; 0; 0];

        % OBSTACLE 3: circle 
        r = 0.3;
        obs_all(:,3)     = [R*(6 - r*cos(w*t)); 0 + r*R*sin(w*t); 0];
        obs_dot_all(:,3) = [r*R*w*sin(w*t);       r*R*w*cos(w*t); 0];
        obs_2dot_all(:,3)= [r*R*(w^2)*cos(w*t);  -r*R*(w^2)*sin(w*t); 0];

end

end
