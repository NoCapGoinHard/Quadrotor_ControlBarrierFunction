%% December 2025
%% Carlo Rugiero, Francesco Maria Germano, Giovanni Pio Cuoco, Matteo Trusiani
%%
%% Control Barrier Function for collision avoidance motion
%% for an unconstrained robot modeled as a double integrator
%% (quadrotor under simplified hierarchical control)

clc; close all; clear variables;

%% parameters
global traj_type M T_s kp kd cbf_type mu delta delta1 V P_obs obs_est obs_dot_est obs_2dot_est obs_3dot_est obs_4dot_est next_print_t

traj_type = 'line';     % 'line' or 'circle' or 'square'

% robot initial conditions. state=[x y z v_x v_y v_z]'
switch traj_type
    case 'line'
        initialConditions=[0;0;0;0;0;0];
        T=20;
    case 'circle'
        initialConditions=[6;0;1;0;0;0];
        T=2*pi;
    case 'square'
        initialConditions=[-4.7;4.7;0.3;0;0;0];
        T=40;
    otherwise
        error('Please select traj_type among the available values');
end

M = 3; % number of obstacles

T_s=0.005; % sampling time

% reference controller parameters
kp = 25; % proportional gain
kd = 10; % derivative gain

% control barrier function (cbf) parameters
cbf_type = 'static'; % 'static' or 'dynamic'
delta1 = 0.2; % cbf activation thereshold
delta = delta1/10; % collision thereshold
mu = 0.5; % cbf gain

% Kalman filter parameters
V=0.001*eye(15*M); % process noise covariance
P_obs=V; % initial value of estimate coviariance
obs_est=zeros(3*M,1); % initial value of obstacle position estimate
obs0 = obs_traj_multi(0);
for i=1:M
    obs_est(3*i-2:3*i)=obs0(:,i);
end
obs_dot_est=zeros(3*M,1); % initial value of obstacle velocity estimate
obs_2dot_est=zeros(3*M,1); % initial value of obstacle acceleration estimate
obs_3dot_est=zeros(3*M,1); % initial value of obstacle acceleration estimate
obs_4dot_est=zeros(3*M,1); % initial value of obstacle acceleration estimate

% printing rate of the debugging string
next_print_t = 0;

% selection_animations tells the function animations() which plots to generate
% selection_animations(1)/=0 plots trajectories (with time)
% selection_animations(2)/=0 plots paths (no time)
% selection_animations(3)/=0 plots minimum distance
selection_animations=[0;10;10];

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
xd   = zeros(length(t),1);
yd   = zeros(length(t),1);
xobs = zeros(length(t),M);
yobs = zeros(length(t),M);
d_min  = zeros(length(t),1);

% compute reference and obstacles trajectories
for i = 1:length(t)
    [pd,~,~,~,~] = traj_plan(t(i));
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

animations(t,x,y,xd,yd,xobs,yobs,d_min,selection_animations);






%% Controller
function u=controller(t,state)
global cbf_type mu delta delta1 kp kd next_print_t M
p = state(1:3); % robot position
p_dot = state(4:6); % robot velocity

% obstacle
obs_all=obs_traj_multi(t);
[obs_dot_all, obs_2dot_all,~,~]=kalman(obs_all);

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
z     = z_all(:,i_star);   % z = p-p_obs,i
switch cbf_type
    case 'dynamic'
        z_dot = z_dot_all(:,i_star);
        obs_2dot = obs_2dot_all(:,i_star);
    case 'static'
        z_dot = p_dot;
        obs_2dot = zeros(3,1);
    otherwise
        error('Please select cbf_type among the available values');
end


% checks for collisions
if norm(z)<=delta
    error(['Collision has happened at t = ', num2str(t)]);
end

% trajectory planning
[pd, pd_dot, pd_2dot,~,~]=traj_plan(t);

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
if t >= next_print_t
    disp(['t = ', num2str(t), ...
        ', h = ', num2str(h), ...
        ', h_dot_star = ', num2str(h_dot_star), ...
        ', coeff = ', num2str(coeff), ...
        '.']);
    next_print_t = next_print_t + 1;
end
end

