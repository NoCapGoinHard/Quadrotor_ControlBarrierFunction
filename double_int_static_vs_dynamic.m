%% December 2025
%% Carlo Rugiero, Francesco Maria Germano, Giovanni Pio Cuoco, Matteo Trusiani
%%
%% Control Barrier Function for collision avoidance motion
%% for an unconstrained robot modeled as a double integrator
%% (quadrotor under simplified hierarchical control)

clc; close all; clear variables;

%% parameters
global obs_type M T_s kp kd cbf_type mu delta delta1 V P_obs obs_est obs_dot_est obs_2dot_est obs_3dot_est obs_4dot_est next_print_t

% robot initial conditions. state=[x y z v_x v_y v_z]'
initialConditions=[0.3;0.3;0;0;0;0];

T_s=0.005; % sampling time
T=10; % total simulation length
M=1;
obs_type = 'parabola'; % 'parabola' or 'point' or 'vert_line' or 'hor_line'

% reference controller parameters
kp = 100; % proportional gain
kd = 20; % derivative gain

% control barrier function (cbf) parameters
cbf_type = 'dynamic'; % 'static' or 'dynamic'
delta1 = 0.15; % cbf activation thereshold
delta = delta1/10; % collision thereshold
mu = 0.5; % cbf gain

% Kalman filter parameters
V=0.001*eye(15); % process noise covariance
P_obs=V; % initial value of estimate coviariance
obs_est=obs_traj(0); % initial value of obstacle position estimate
obs_dot_est=zeros(3,1); % initial value of obstacle velocity estimate
obs_2dot_est=zeros(3,1); % initial value of obstacle acceleration estimate
obs_3dot_est=zeros(3,1); % initial value of obstacle acceleration estimate
obs_4dot_est=zeros(3,1); % initial value of obstacle acceleration estimate

% printing rate of the debugging string
next_print_t = 0;

% selection_animations tells the function animations() which plots to generate
% selection_animations(1)/=0 plots trajectories (with time)
% selection_animations(2)/=0 plots paths (no time)
% selection_animations(3)/=0 plots minimum distance
selection_animations=[10;0;0];

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
xobs = zeros(length(t),1);
yobs = zeros(length(t),1);

% compute reference and obstacles trajectories
for i = 1:length(t)
    xd(i) = t(i);
    yd(i) = 0;
    obs = obs_traj(t(i));
    xobs(i)=obs(1);
    yobs(i)=obs(2);
end

animations(t,x,y,xd,yd,xobs,yobs,selection_animations);





%% Controller
function u=controller(t,state)
global cbf_type mu delta delta1 kp kd next_print_t M
p = state(1:3); % robot position
p_dot = state(4:6); % robot velocity

% obstacle
obs=obs_traj(t);
[obs_dot, obs_2dot,~,~]=kalman(obs);
z = p-obs;
switch cbf_type
    case 'dynamic'
        z_dot = p_dot-obs_dot;
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
 pd      = [t; 0; 0]; 
 pd_dot  = [1; 0; 0]; 
 pd_2dot = [0; 0; 0];

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

%% Obstacle
function obs=obs_traj(t)
global obs_type

switch obs_type
    case 'hor_line'
        vel=1;
        obs=[-vel*t+5*(1+vel); 0; 0];
    case 'vert_line'
        vel=1;
        obs=[5; vel*(-t+5); 0];
    case 'point'
        obs=[5; 0; 0];
    case 'parabola'
        obs=[t; (t-3)*(t-7); 0];
    otherwise
        error('Please select obs_type among the available values');
end
end