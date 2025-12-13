%% December 2025
%% Carlo Rugiero, Francesco Maria Germano, Giovanni Pio Cuoco, Matteo Trusiani
%%
%% Control Barrier Function for collision avoidance motion
%% for an unconstrained robot modeled as a double integrator
%% (quadrotor under Dynamic Feedback Linearization)

clc; close all; clear variables;

%% parameters
global traj_type M T_s R w k0 k1 k2 k3 mu2 mu3 mu4 delta delta1 V P_obs obs_est obs_dot_est obs_2dot_est obs_3dot_est obs_4dot_est next_print_t

% robot initial conditions. state=[x y z  v_x v_y v_z  a_x a_y a_z  j_x j_y j_z]'
initialConditions=[0.3;0.3;0; 0;0;0; 0;0;0; 0;0;0];

traj_type = 'line';     % 'line' or 'circle' or 'square'

M = 3; % number of obstacles

R = 1; % radious of the circumpherence
w = 1; % angular velocity

T_s=0.005; % sampling time

% reference controller parameters
if strcmp(traj_type,'circle')
    eig=1;
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

% selection_vector tells the function animations() which plots to generate
% selection_vector(1)/=0 plots trajectories (with time)
% selection_vector(2)/=0 plots paths (no time)
% selection_vector(3)/=0 plots minimum distance
selection_vector=[0;10;0];

%% running ode

% total simulation length
switch traj_type
    case 'line'
        T=20;
    case 'circle'
        T = 2*pi/w;
    case 'square'
        T = 40;
    otherwise
        error('Please select traj_type among the available values');
end

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

animations(t,x,y,xd,yd,xobs,yobs,selection_vector);





%% Controller
function u=controller(t,state)
global mu2 mu3 mu4 delta delta1 k0 k1 k2 k3 next_print_t M
p = state(1:3); % robot position
p_dot = state(4:6); % robot velocity
p_2dot = state(7:9); % robot velocity
p_3dot = state(10:12); % robot velocity

% obstacles
obs_all=obs_traj_multi(t);
[obs_dot_all, obs_2dot_all, obs_3dot_all, obs_4dot_all]=kalman(obs_all);

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
if t >= next_print_t
    disp(['t = ', num2str(t), ...
        ', h = ', num2str(h), ...
        ', h_dot_star = ', num2str(h_dot_star), ...
        ', coeff = ', num2str(coeff), ...
        '.']);
    next_print_t = next_print_t + 1;
end
end

