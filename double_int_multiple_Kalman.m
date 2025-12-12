%% November 2025
%% Carlo Rugiero, Francesco Maria Germano, Giovanni Pio Cuoco, Matteo Trusiani

%% Control Barrier Function for collision avoidance motion
%% for an unconstrained robot modeled as a double integrator
%% (quadrotor under simplified hierarchical control)

clc; close all; clear variables;

%% parameters
global traj_type M T_s R w kp kd mu delta delta1 V P_obs obs_est obs_dot_est obs_2dot_est next_print_t_2

% robot initial conditions. state=[x y z v_x v_y v_z]'
initialConditions=[0.3;0.3;0;0;0;0];

traj_type = 'square';     % 'line' or 'circle' or 'square'

M = 3; % number of obstacles

% check condition to avoid overlapping
obs0 = obs_traj_multi(0);
for i = 1:M
    for j = i+1:M
        if norm(obs0(:,i) - obs0(:,j))^2 <= 2*delta1
            error('Obstacles violate condition ||p_i - p_j||^2 > 2 delta1');
        end
    end
end

R = 1; % radious of the circumpherence
w = 1; % angular velocity

% T_s: sampling time, T: total simulation length
T_s=0.005;

switch traj_type
    case 'line'
        T=20;
    case 'circle'
        T = 2*pi/w;
    case 'square'
        T = 40;
end

% reference controller parameters
kp = 100; % proportional gain
kd = 20; % derivative gain

% control barrier function (cbf) parameters
delta1 = 0.15; % cbf activation thereshold
delta = delta1/10; % collision thereshold
mu = 0.5; % cbf gain

% Kalman filter parameters
V=0.001*eye(9*M); % process noise covariance
P_obs=V; % initial value of estimate coviariance
obs_est=zeros(3*M,1); % initial value of obstacle position estimate
for i=1:M
    obs_est(3*i-2:3*i)=obs0(:,i);
end
obs_dot_est=zeros(3*M,1); % initial value of obstacle velocity estimate
obs_2dot_est=zeros(3*M,1); % initial value of obstacle acceleration estimate

% printing rate of the debugging string
next_print_t_2 = 1;

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
animations(t,state,M);

%% get the results
x = state(:,1);
y = state(:,2);
xd   = zeros(length(t),1);
yd   = zeros(length(t),1);
xobs = zeros(length(t),M);
yobs = zeros(length(t),M);
d_min  = zeros(length(t),1);

% compute reference and obstacles trajectories
for i = 1:length(t)
    [pd, ~, ~] = traj_plan(t(i));
    xd(i) = pd(1);
    yd(i) = pd(2);
    obs_all = obs_traj_multi(t(i));
    dist2 = zeros(1,M);
    for j = 1:M
        xobs(i,j) = obs_all(1,j); % the x-coordinate of obstacle j at time t(i)
        yobs(i,j) = obs_all(2,j); % the y-coordinate of obstacle j at time t(i)
        dist2(j)  = (state(i,1) - obs_all(1,j))^2 + (state(i,2) - obs_all(2,j))^2; % squared Euclidean distance = (x_robot - x_obstacle)^2 + (y_robot - y_obstacle)^2
    end
    % minimum distance
    d_min(i) = sqrt(min(dist2));
end

% plot robot and obstacles paths (no time)
% figure(1); hold on; grid on; axis equal;
% plot(x, y, 'Color', [0 0 0.55], 'LineWidth', 1.8);
% plot(xd, yd, '--', 'Color', [1 0.4 0.2], 'LineWidth', 1.5);
% for j = 1:M
%     plot(xobs(:,j), yobs(:,j), '--','Color',[0 0 0], 'LineWidth', 2);
% end
% legend('Robot', 'Reference', 'Obstacles');
% xlabel('x [m]'); ylabel('y [m]');
% title('Robot trajectory with multiple obstacles');

% plot robot and obstacles trajectories (with time)
switch traj_type
    case 'line'
        x_limit = [-1 21];
        y_limit = [-1 21];
        step_size = 1000;
    case 'circle'
        x_limit = [-2 2];
        y_limit = [-2 2];
        step_size = 100;
    case 'square'
       x_limit = [-5 15];
       y_limit = [-15 5];
       step_size = 1000;
end
num_loops = 2;
tail_length = 100;
N = length(t);
indices = 1:step_size:N;

figure(1); hold on; grid on; axis equal;
xlabel('x [m]'); ylabel('y [m]');
title('Robot Trajectory with Multiple Dynamic Obstacles');
axis([x_limit(1), x_limit(2), y_limit(1), y_limit(2)]); 
h_ref_line = plot(xd, yd, '--', 'Color', [1 0.4 0.2], 'LineWidth', 1.5);
h_desired_marker = plot(xd(1), yd(1), 'd', 'MarkerSize', 8, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
h_robot_path = plot(x(1), y(1), 'Color', [0 0 0.55], 'LineWidth', 1.8);
h_obs_marker = zeros(M, 1);
for j = 1:M
h_obs_marker(j) = plot(xobs(1,j), yobs(1,j), 's', 'MarkerSize', 8, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
end
h_robot_marker = plot(x(1), y(1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0 0 0.55], 'MarkerEdgeColor', [0 0 0.55]);
h_time = text(x(1), y(1) + 0.2, ['t = ' num2str(t(1), '%.2f') ' s'], 'FontSize', 10);
legend([h_robot_marker, h_robot_path, h_ref_line, h_desired_marker, h_obs_marker(1)], ...
'Robot Position', 'Robot Path', 'Reference Line', 'Desired Position', 'Obstacles', ...
'Location', 'best');
for loop_count = 1:num_loops
t_prev = t(1);
set(h_robot_path, 'XData', x(1), 'YData', y(1));
for idx = indices(2:end)
start_idx = max(1, idx - tail_length); 
dt = t(idx) - t_prev;
set(h_robot_path, 'XData', x(start_idx:idx), 'YData', y(start_idx:idx));
set(h_desired_marker, 'XData', xd(idx), 'YData', yd(idx));
set(h_robot_marker, 'XData', x(idx), 'YData', y(idx));
for j = 1:M
set(h_obs_marker(j), 'XData', xobs(idx, j), 'YData', yobs(idx, j));
end
set(h_time, 'String', ['t = ' num2str(t(idx), '%.2f') ' s']);
set(h_time, 'Position', [x(idx), y(idx) + 0.2]);
drawnow;
pause(dt);
t_prev = t(idx);
end
pause(0.5); 
end
hold off;

% plot minimum distance vs time
% figure(2); hold on; grid on;
% plot(t, d_min, 'b', 'LineWidth', 1.8);
% yline(delta,  'r--', 'LineWidth', 2);
% yline(delta1, 'g--', 'LineWidth', 2);
% xlabel('time [s]');
% ylabel('min distance [m]');
% legend('minimum distance', 'collision threshold \delta', 'CBF threshold \delta_1');
% title('Minimum distance to obstacles over time');








%% Controller
function u=controller(t,state)
global mu delta delta1 kp kd next_print_t_2 M
p = state(1:3); % robot position
p_dot = state(4:6); % robot velocity

% obstacle
obs_all=obs_traj_multi(t);
[obs_dot_all, obs_2dot_all]=kalman(obs_all);

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
obs_2dot = obs_2dot_all(:,i_star);
z     = z_all(:,i_star);   % z = p-p_obs,i
z_dot = z_dot_all(:,i_star);

% checks for collisions
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
    next_print_t_2 = next_print_t_2 + 1;
end
end

%% Trajectory planner
function [pd, pd_dot, pd_2dot]=traj_plan(t)
% generates a reference trajectory to follow
global traj_type R w
switch traj_type
    case 'line'
        pd      = [t; 0; 0]; pd_dot  = [1; 0; 0]; pd_2dot = [0; 0; 0];
    case 'circle'
        pd      = [       R*cos(w*t);        R*sin(w*t); 0];
        pd_dot  = [    -R*w*sin(w*t);      R*w*cos(w*t); 0];
        pd_2dot = [-R*(w^2)*cos(w*t); -R*(w^2)*sin(w*t); 0];
    case 'square'
        if t<=10
            pd      = [t; 0; 0]; pd_dot  = [1; 0; 0]; pd_2dot = [0; 0; 0];
        elseif t<=20
            pd      = [10; 10-t; 0]; pd_dot  = [0; -1; 0]; pd_2dot = [0; 0; 0];
        elseif t<=30
            pd      = [30-t; -10; 0]; pd_dot  = [-1; 0; 0]; pd_2dot = [0; 0; 0];
        elseif t<=40
            pd      = [0; -40+t; 0]; pd_dot  = [0; 1; 0]; pd_2dot = [0; 0; 0];
        end
    otherwise
        error('Please select ref_type among the available values');
end
end

%% Obstacles
function obs_all=obs_traj_multi(t)
% generates the obstacle trajectory
global traj_type w R M
obs_all     = zeros(3,M);
switch traj_type
    case 'line'
        % OBSTACLE 1: line
        vel=5;
        alpha=pi/10;
        obs_all(:,1)=[-vel*t+5*(1+vel); vel*(t-5)*tan(alpha); 0];
        % OBSTACLE 2: static
        obs_all(:,2)     = [9; 0; 0];
        % OBSTACLE 3: parabola
        obs_all(:,2)=[t; t^2-27*t+180; 0];
    case 'circle'
        % OBSTACLE 1: static
        obs_all(:,1)     = [0; -1; 0];
        % OBSTACLE 2: vertical line
        obs_all(:,2)     = [-1; pi - t; 0];
        % OBSTACLE 3: circle
        r = 0.3;
        obs_all(:,3)     = [R*(0 - r*cos(w*t)); 1 + r*R*sin(w*t); 0];
    case 'square'
        % OBSTACLE 1: circle
        r = 3;
        obs_all(:,1)     = [r*cos(-w*t+3.1); r*sin(-w*t+3.1); 0];
        % OBSTACLE 2: parabola
        obs_all(:,2)     = [t^2-30*t+231; 10-t; 0];
        % OBSTACLE 3: line
        obs_all(:,3)     = [-1+5*cos(t); -6-5*cos(t); 0];
end
end

%% Kalman filter
function [obs_dot, obs_2dot]=kalman(obs)
global T_s P_obs V obs_est obs_dot_est obs_2dot_est M

% discrete model of the obstacle motion
A = [eye(3*M), T_s*eye(3*M), ((T_s^2)/2)*eye(3*M);
    zeros(3*M), eye(3*M), T_s*eye(3*M);
    zeros(3*M), zeros(3*M), eye(3*M)];
C = [eye(3*M) zeros([3*M 6*M])];

W=0.01*eye(3*M);  % measurement noise

% prediction step
x_pred=A*[obs_est;obs_dot_est;obs_2dot_est];
P_pred=A*P_obs*A'+V;
Gain=P_pred*C'/(C*P_pred*C'+W); % Kalman gain

% correction step
obs_vector=zeros(3*M,1);
for i=1:M
    obs_vector(3*i-2:3*i)=obs(:,i);
end
Inn=obs_vector-C*x_pred; % innovation
x_corr=x_pred+Gain*Inn;
P_obs=(eye(9*M)-Gain*C)*P_pred;

% update variables for the next iteration
obs_est=x_corr(1:3*M);
obs_dot_est=x_corr(3*M+1:6*M);
obs_2dot_est=x_corr(6*M+1:9*M);

obs_dot=zeros(3,M);
obs_2dot=zeros(3,M);
% output of the function
for i=1:M
    obs_dot(:,i)=obs_dot_est(3*i-2:3*i,1);
    obs_2dot(:,i)=obs_2dot_est(3*i-2:3*i,1);
end

end