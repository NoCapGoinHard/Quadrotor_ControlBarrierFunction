%% November 2025
%% Carlo Rugiero, Francesco Maria Germano, Giovanni Pio Cuoco, Matteo Trusiani

%% Control Barrier Function for collision avoidance motion
%% for an unconstrained robot modeled as a double integrator
%% (quadrotor under Dynamic Feedback Linearization)

clc; close all; clear variables;

%% parameters
global traj_type M T_s R w k0 k1 k2 k3 mu2 mu3 mu4 delta delta1 V P_obs obs_est obs_dot_est obs_2dot_est obs_3dot_est obs_4dot_est next_print_t_1 next_print_t_2

% robot initial conditions. state=[x y z  v_x v_y v_z  a_x a_y a_z  j_x j_y j_z]'
initialConditions=[1.25;0.25;0; 0;0;0; 0;0;0; 0;0;0];

traj_type = 'line';     % 'line' or 'circle'

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

% T_s: sampling time, T: total simulation length
T_s=0.005; 
if strcmp(traj_type,'circle')
    R = 1; % radious of the circumpherence
    w = 1; % angular velocity
    T = 2*pi/w;
else
    T = 15;
end

% reference controller parameters
if strcmp(traj_type,'circle')
    eig=3.5;
else
    eig=10;
end
k0 = 1*eig^4;
k1 = 4*eig^3;
k2 = 6*eig^2;
k3 = 4*eig;

% control barrier function (cbf) parameters
delta1 = 0.1; % cbf activation thereshold
delta = delta1/10; % collision thereshold
mu=1;
q=10;
mu2 = mu*q; 
mu3 = 0;
mu4 = mu;

% Kalman filter parameters
V=0.001*eye(15*M); % process noise covariance
P_obs=V; % initial value of estimate coviariance
obs_est=zeros(3*M,1); % initial value of obstacle position estimate
for i=1:M
obs_est(3*i-2:3*i)=obs0(:,i); 
end
obs_dot_est=zeros(3*M,1); % initial value of obstacle velocity estimate
obs_2dot_est=zeros(3*M,1); % initial value of obstacle acceleration estimate
obs_3dot_est=zeros(3*M,1); % initial value of obstacle acceleration estimate
obs_4dot_est=zeros(3*M,1); % initial value of obstacle acceleration estimate

% printing rate of the debugging string
next_print_t_1 = 0.1;
next_print_t_2 = 0.1;

%% running ode
state_current = initialConditions;
t_current = 0;
steps = 1:T/T_s;
t=[]; state=[];

for i=steps
    u=controller(t_current,state_current);
    [t_ode,state_ode]=ode45(@(t_ode,state_ode) [state_ode(4:6);state_ode(7:9);state_ode(10:12);u], [(i-1)*T_s (i)*T_s], state_current, odeset('RelTol',1e-9,'AbsTol',1e-15));

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

% plot robot and obstacles trajectories
figure(1); hold on; grid on; axis equal;
plot(x, y, 'Color', [0 0 0.55], 'LineWidth', 1.8);
plot(xd, yd, '--', 'Color', [1 0.4 0.2], 'LineWidth', 1.5);
for j = 1:M
    plot(xobs(:,j), yobs(:,j), '--','Color',[0 0 0], 'LineWidth', 2);
end
legend('Robot', 'Reference', 'Obstacles');
xlabel('x [m]'); ylabel('y [m]');
title('Robot trajectory with multiple obstacles');

% plot minimum distance vs time
figure(2); hold on; grid on;
plot(t, d_min, 'b', 'LineWidth', 1.8);
yline(delta,  'r--', 'LineWidth', 2);
yline(delta1, 'g--', 'LineWidth', 2);
xlabel('time [s]');
ylabel('min distance [m]');
legend('minimum distance', 'collision threshold \delta', 'CBF threshold \delta_1');
title('Minimum distance to obstacles over time');

%% Controller
function u=controller(t,state)
global mu2 mu3 mu4 delta delta1 k0 k1 k2 k3 next_print_t_2 M
p = state(1:3); % robot position
p_dot = state(4:6); % robot velocity
p_2dot = state(7:9); % robot velocity
p_3dot = state(10:12); % robot velocity

% obstacle
obs_all=obs_traj_multi(t);
[obs_dot_all, obs_2dot_all, obs_3dot_all, obs_4dot_all]=kalman(t,obs_all);

% compute relative vectors and distances
z_all     = zeros(3,M);
z_dot_all = zeros(3,M);
z_2dot_all = zeros(3,M);
z_3dot_all = zeros(3,M);
dist2     = zeros(1,M);

for i = 1:M
    z_all(:,i)     = p     - obs_all(:,i);
    z_dot_all(:,i) = p_dot - obs_dot_all(:,i);
    z_2dot_all(:,i) = p_2dot - obs_2dot_all(:,i);
    z_3dot_all(:,i) = p_3dot - obs_3dot_all(:,i);
    dist2(i)       = z_all(:,i)'*z_all(:,i); % = ||p - p_obs,i||^2
end

% select the closest obstacle
[~, i_star] = min(dist2);
obs_4dot = obs_4dot_all(:,i_star);
z     = z_all(:,i_star);   % z = p-p_obs,i
z_dot = z_dot_all(:,i_star);
z_2dot = z_2dot_all(:,i_star);
z_3dot = z_3dot_all(:,i_star);

% checks for collisions
if norm(z)<=delta
    error(['Collision has happened at t = ', num2str(t)]);
end

% trajectory planning
[pd, pd_dot, pd_2dot, pd_3dot, pd_4dot]=traj_plan(t);

% trajectory tracking controller (reference controller in absence of obstacles)
u_star = pd_4dot + k3*(pd_3dot-p_3dot) + k2*(pd_2dot-p_2dot) + k1*(pd_dot-p_dot) + k0*(pd-p);

% control barrier function (controller with obstacles)
h=z'*(z+mu2*z_dot+mu3*z_2dot+mu4*z_3dot);
h_dot_star = z_dot'*(z+mu2*z_dot+mu3*z_2dot+mu4*z_3dot)+z'*(z_dot+mu2*z_2dot+mu3*z_3dot+mu4*(u_star-obs_4dot));

if (h<=delta1)&&(h_dot_star<=0)
    u = obs_4dot -(2*z_dot+mu2*z_2dot+mu3*z_3dot)/mu4 -(mu3*z_dot'*z_2dot+mu4*z_dot'*z_3dot)*(z+mu2*z_dot+mu3*z_2dot+mu4*z_3dot)/(mu4*h) +(eye(3)-(z*z')/(z'*z))*u_star;
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

%% Trajectory planner
function [pd, pd_dot, pd_2dot, pd_3dot, pd_4dot]=traj_plan(t)
% generates a reference trajectory to follow
global traj_type R w
switch traj_type
    case 'line'
        pd      = [t; 0; 0];
        pd_dot  = [1; 0; 0];
        pd_2dot = [0; 0; 0];
        pd_3dot = [0; 0; 0];
        pd_4dot = [0; 0; 0];
    case 'circle'
        pd      = [       R*cos(w*t);        R*sin(w*t); 0];
        pd_dot  = [    -R*w*sin(w*t);      R*w*cos(w*t); 0];
        pd_2dot = [-R*(w^2)*cos(w*t); -R*(w^2)*sin(w*t); 0];
        pd_3dot = [ R*(w^3)*sin(w*t); -R*(w^3)*cos(w*t); 0];
        pd_4dot = [ R*(w^4)*cos(w*t);  R*(w^4)*sin(w*t); 0];
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
    case 'circle'
        % OBSTACLE 1: static
        obs_all(:,1)     = [0; -1; 0];
        % OBSTACLE 2: vertical line
        obs_all(:,2)     = [-1; pi - t; 0];
        % OBSTACLE 3: circle
        r = 0.3;
        obs_all(:,3)     = [R*(0 - r*cos(w*t)); 1 + r*R*sin(w*t); 0];
    case 'line'
        % OBSTACLE 1: static
        obs_all(:,1)     = [4; 0; 0];
        % OBSTACLE 2: vertical line
        obs_all(:,2)     = [3; pi - t; 0];
        % OBSTACLE 3: circle
        r = 0.3;
        obs_all(:,3)     = [R*(6 - r*cos(w*t)); 0 + r*R*sin(w*t); 0];
end
end

%% Kalman filter
function [obs_dot, obs_2dot, obs_3dot, obs_4dot]=kalman(t,obs)
global T_s P_obs V obs_est obs_dot_est obs_2dot_est obs_3dot_est obs_4dot_est next_print_t_1 M

% discrete model of the obstacle motion
A = [eye(3*M), T_s*eye(3*M), ((T_s^2)/2)*eye(3*M), ((T_s^3)/6)*eye(3*M), ((T_s^4)/24)*eye(3*M);
    zeros(3*M), eye(3*M), T_s*eye(3*M), ((T_s^2)/2)*eye(3*M), ((T_s^3)/6)*eye(3*M);
    zeros(3*M), zeros(3*M), eye(3*M), T_s*eye(3*M), ((T_s^2)/2)*eye(3*M);
    zeros(3*M), zeros(3*M), zeros(3*M), eye(3*M), T_s*eye(3*M);
    zeros(3*M), zeros(3*M), zeros(3*M), zeros(3*M), eye(3*M)];
C = [eye(3*M) zeros([3*M 12*M])];

W=0.01*eye(3*M);  % measurement noise

% prediction step
x_pred=A*[obs_est;obs_dot_est;obs_2dot_est;obs_3dot_est;obs_4dot_est];
P_pred=A*P_obs*A'+V;
Gain=P_pred*C'/(C*P_pred*C'+W); % Kalman gain

% correction step
obs_vector=zeros(3*M,1);
for i=1:M
obs_vector(3*i-2:3*i)=obs(:,i); 
end
Inn=obs_vector-C*x_pred; % innovation
x_corr=x_pred+Gain*Inn;
P_obs=(eye(15*M)-Gain*C)*P_pred;

% update variables for the next iteration
obs_est=x_corr(1:3*M);
obs_dot_est=x_corr(3*M+1:6*M);
obs_2dot_est=x_corr(6*M+1:9*M);
obs_3dot_est=x_corr(9*M+1:12*M);
obs_4dot_est=x_corr(12*M+1:15*M);

% print prediction error for debugging
if t >= next_print_t_1
    disp(['|est_error| = ',num2str(norm(obs_vector-obs_est))]);
    next_print_t_1 = next_print_t_1 + 0.1;
end

obs_dot=zeros(3,M);
obs_2dot=zeros(3,M);
obs_3dot=zeros(3,M);
obs_4dot=zeros(3,M);
% output of the function
for i=1:M
obs_dot(:,i)=obs_dot_est(3*i-2:3*i,1);
obs_2dot(:,i)=obs_2dot_est(3*i-2:3*i,1);
obs_3dot(:,i)=obs_3dot_est(3*i-2:3*i,1);
obs_4dot(:,i)=obs_4dot_est(3*i-2:3*i,1);
end

end