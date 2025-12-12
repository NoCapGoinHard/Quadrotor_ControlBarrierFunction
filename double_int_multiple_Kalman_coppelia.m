%% November 2025
%% Carlo Rugiero, Francesco Maria Germano, Giovanni Pio Cuoco, Matteo Trusiani

%% Control Barrier Function for collision avoidance motion
%% for an unconstrained robot modeled as a double integrator
%% (quadrotor under simplified hierarchical control)

clc; close all; clear variables;

%% Connect to CoppeliaSim (ZeroMQ Remote API)
client = RemoteAPIClient();
sim = client.getObject('sim');
%% DEBUG COPPELIASIM INTERACTION
fprintf('\n[DBG] Connected. sim class: %s\n', class(sim));

% Basic sanity checks on sim object
dbgPrintSimState(sim, 'BEFORE startSimulation');

% Ensure stepping
try
    sim.setStepping(true);
    fprintf('[DBG] setStepping(true) OK\n');
catch ME
    fprintf('[DBG][ERR] setStepping(true) failed: %s\n', ME.message);
end

% Start simulation
try
    sim.startSimulation();
    fprintf('[DBG] startSimulation() OK\n');
catch ME
    fprintf('[DBG][ERR] startSimulation() failed: %s\n', ME.message);
end

% Do a few warm-up steps to let sysCall_init/sysCall_sensing publish signals
for k = 1:5
    try
        sim.step();
        fprintf('[DBG] step %d OK\n', k);
    catch ME
        fprintf('[DBG][ERR] step %d failed: %s\n', k, ME.message);
    end
end

dbgPrintSimState(sim, 'AFTER startSimulation + warmup');

% Optional: write test inputs (so Lua reads something non-nil)
safeSetFloatSignal(sim, 'ControlInputs/ax', 0.0);
safeSetFloatSignal(sim, 'ControlInputs/ay', 0.0);
safeSetFloatSignal(sim, 'ControlInputs/az', 0.0);

% Now attempt to read the state signals with safe wrapper
x_test  = safeGetFloatSignal(sim, 'State/x');
y_test  = safeGetFloatSignal(sim, 'State/y');
z_test  = safeGetFloatSignal(sim, 'State/z');
vx_test = safeGetFloatSignal(sim, 'State/vx');
vy_test = safeGetFloatSignal(sim, 'State/vy');
vz_test = safeGetFloatSignal(sim, 'State/vz');

fprintf('[DBG] Initial read summary: x=%.6f y=%.6f z=%.6f vx=%.6f vy=%.6f vz=%.6f\n', ...
    x_test, y_test, z_test, vx_test, vy_test, vz_test);

% Try listing float signals (may be unsupported)
listFloatSignalsIfPossible(sim);

% If still NaN, stop immediately (do not continue into control loop)
if any(isnan([x_test y_test z_test vx_test vy_test vz_test]))
    error('[DBG] State signals still not readable. Stopping here for diagnosis.');
end

%% parameters
global traj_type M T_s R w kp kd mu delta delta1 V P_obs obs_est obs_dot_est obs_2dot_est next_print_t_2

% robot initial conditions. state=[x y z v_x v_y v_z]'
%initialConditions=[0;0;0;0;0;0]; %NOW SET UP IN LUA SCRIPT

traj_type = 'line';     % 'line' or 'circle' or 'square'

M = 3; % number of obstacles

% control barrier function (cbf) parameters
delta1 = 0.15; % cbf activation thereshold
delta = delta1/10; % collision thereshold
mu = 0.5; % cbf gain

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
dt_scene = sim.getSimulationTimeStep();
T_s= dt_scene; 

if strcmp(traj_type,'circle')
    T = 2*pi/w;
elseif strcmp(traj_type,'square')
    T = 40;
else
    T=20;
end

% reference controller parameters
kp = 100; % proportional gain
kd = 20; % derivative gain



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


%% WHILE LOOP INSTEAD OF ODE (for CoppeliaSim)
% real-time control loop
sim.setStepping(true);
sim.startSimulation();

% --- Wait until State signals are available (handshake) ---
maxWaitSteps = 200; % 200 steps * dt (scene) e.g. 200*0.05 = 10s max
ok = false;

for k = 1:maxWaitSteps
    sim.step();
    x_try = safeGetFloatSignal(sim, 'State/x'); % your safe helper
    if ~isnan(x_try)
        fprintf('[DBG] State/x available after %d steps\n', k);
        ok = true;
        break;
    end
end

if ~ok
    error('[DBG] State signals not published by Lua. Check script type/errors and signal names.');
end


t = 0;
log_t = [];
log_state = [];

while t < T
    sim.step();
    % ---- READ STATE FROM COPPELIASIM ----
    x  = sim.getFloatSignal('State/x');
    y  = sim.getFloatSignal('State/y');
    z  = sim.getFloatSignal('State/z');
    
    vx = sim.getFloatSignal('State/vx');
    vy = sim.getFloatSignal('State/vy');
    vz = sim.getFloatSignal('State/vz');
    
    state = [x; y; z; vx; vy; vz];
    
    u = controller(t, state);
    
    sim.setFloatSignal('ControlInputs/ax', u(1));
    sim.setFloatSignal('ControlInputs/ay', u(2));
    sim.setFloatSignal('ControlInputs/az', u(3));

    log_t(end+1,1) = t;
    log_state(end+1,:) = state.';

    % ---- TIME UPDATE ----
    t = t + T_s;

end

sim.stopSimulation();

%% get the results
t_vec = log_t;

x = log_state(:,1);
y = log_state(:,2);

xd   = zeros(length(t_vec),1);
yd   = zeros(length(t_vec),1);
xobs = zeros(length(t_vec),M);
yobs = zeros(length(t_vec),M);
d_min  = zeros(length(t_vec),1);

for i = 1:length(t_vec)

    ti = t_vec(i);

    % reference
    [pd, ~, ~] = traj_plan(ti);
    xd(i) = pd(1);
    yd(i) = pd(2);

    % obstacles
    obs_all = obs_traj_multi(ti);

    dist2 = zeros(1,M);
    for j = 1:M
        xobs(i,j) = obs_all(1,j);
        yobs(i,j) = obs_all(2,j);
        dist2(j)  = (x(i)-obs_all(1,j))^2 + (y(i)-obs_all(2,j))^2;
    end

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

% % print the current relevant datas for debugging
% if t >= next_print_t_2
%     disp(['t = ', num2str(t), ...
%         ', h = ', num2str(h), ...
%         ', h_dot_star = ', num2str(h_dot_star), ...
%         ', coeff = ', num2str(coeff), ...
%         '.']);
%     next_print_t_2 = next_print_t_2 + 1;
% end

if mod(round(t/T_s),200)==0
    disp(['t=',num2str(t),'  u=',num2str(u.')]);
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
         % obs_all(:,3)     = [R*(6 - r*cos(w*t)); 0 + r*R*sin(w*t); 0];
        vel=5;
        alpha=pi/10;
        obs_all(:,3)=[-vel*t+5*(1+vel); vel*(t-5)*tan(alpha); 0];
        % OBSTACLE 1: static
        obs_all(:,1)     = [9; 0; 0];
        % OBSTACLE 2: vertical line
        %obs_all(:,2)     = [3; pi - t; 0];
        obs_all(:,2)=[t; t^2-27*t+180; 0];
        % OBSTACLE 3: circle
        % r = 0.3;
       
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

%% HELPER FUNCTIONS FOR COPPELIA DEBUGGING
function v = safeGetFloatSignal(sim, name)
%SAFEGETFLOATSIGNAL Read a float signal with detailed debug info.
% Returns NaN if it cannot be read.

try
    v = sim.getFloatSignal(name);
    if isempty(v)
        fprintf('[DBG] getFloatSignal("%s") returned EMPTY\n', name);
        v = NaN;
    else
        fprintf('[DBG] getFloatSignal("%s") = %.6f\n', name, v);
    end
catch ME
    fprintf('[DBG][ERR] getFloatSignal("%s") failed: %s\n', name, ME.message);
    v = NaN;
end
end

function ok = safeSetFloatSignal(sim, name, value)
%SAFESETFLOATSIGNAL Set a float signal with debug prints.
ok = true;
try
    sim.setFloatSignal(name, value);
    fprintf('[DBG] setFloatSignal("%s") = %.6f\n', name, value);
catch ME
    fprintf('[DBG][ERR] setFloatSignal("%s") failed: %s\n', name, ME.message);
    ok = false;
end
end

function dbgPrintSimState(sim, tag)
%DBGPRINTSIMSTATE Print simulation state if available.
fprintf('\n========== [DBG] %s ==========\n', tag);
try
    st = sim.getSimulationState();
    fprintf('[DBG] sim.getSimulationState() = %d\n', st);
catch ME
    fprintf('[DBG][WARN] sim.getSimulationState() not available: %s\n', ME.message);
end

try
    t = sim.getSimulationTime();
    fprintf('[DBG] sim.getSimulationTime() = %.6f\n', t);
catch ME
    fprintf('[DBG][WARN] sim.getSimulationTime() not available: %s\n', ME.message);
end

try
    dt = sim.getSimulationTimeStep();
    fprintf('[DBG] sim.getSimulationTimeStep() = %.6f\n', dt);
catch ME
    fprintf('[DBG][WARN] sim.getSimulationTimeStep() not available: %s\n', ME.message);
end
fprintf('================================\n\n');
end

function listFloatSignalsIfPossible(sim)
%LISTFLOATSIGNALSIFPOSSIBLE Try to list float signals (depends on CoppeliaSim version).
try
    sigs = sim.getStringSignal('__nonexistent__'); %#ok<NASGU>
catch
    % ignore
end

try
    names = sim.getFloatSignalNameList(); % This may NOT exist in your version
    fprintf('[DBG] Float signal names (%d):\n', numel(names));
    for i = 1:numel(names)
        fprintf('  - %s\n', names{i});
    end
catch ME
    fprintf('[DBG][WARN] Cannot list float signal names (API may not support it): %s\n', ME.message);
end
end
