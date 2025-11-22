%% Dynamic Feedback Linearization
%% of the quadrotor complete model
%% software developed by Paolo Forni for his bachelor degree thesis
%% supervisor: Marilena Vendittelli
%% June 11, 2012

clear variables; close all;

%% Parameters

% quadrotor
global Ix Iy Iz g m J o c0 c1 c2 c3 REF_CASE
T = 1.8;          %% (planning) period T
g = 9.81;
m = 4.34;
Ix = 0.0820;
Iy = 0.0845;
Iz = 0.1377;
o = 2*pi/T;

J = [Ix 0 0; 0 Iy 0;0 0 Iz];

% Hurwitz polynomial coefficients for controller
c3 = 26;c2 = 253;c1 = 1092; c0 = 1764;

%Initial condition on integrators
% 1 = trivial, 2 = slightly different, 3 = not so much aggressive, 4 =
% aggressive
IC_CASE = 1;

%Reference trajectory case:
% 1 = hover, 2 = circumference, 3 piecewise polynomial
REF_CASE = 1;
%% Planning (simulation length from 0 to n seconds)

Tspan = [0 21];

%% Initial condition on integrators

switch IC_CASE
    case 1
        %%%%%%%%% trivial initial conditions %%%%%%
        x_initial = [0 0 0]';
        v_initial = [0 0 0]';
        roll_initial = 0;
        pitch_initial = 0;
        yaw_initial = 0;
        omega_initial = [0 0 0]';
        f_initial = 0.1;
        f_dot_initial = 0;
    
    case 2
        %%%%%%%%% slightly different initial conditions %%%%%%
        x_initial = [0.2 0.2 -0.2]';
        v_initial = [0 0 0]';
        roll_initial = pi/8;
        pitch_initial = pi/8;
        yaw_initial = -pi/8;
        omega_initial = [0 0 0]';
        f_initial = 0.1;
        f_dot_initial = 0;
    case 3
        % %%%%%%%% not-so-much-aggressive initial conditions %%%%
        x_initial = [1 1 2]';
        v_initial = [0 0 0]';
        roll_initial = pi/4;
        pitch_initial = pi/4;
        yaw_initial = pi/4;
        omega_initial = [0 0 0]';
        f_initial = 0.1;
        f_dot_initial = 0;
    case 4
        % %%%%%%%% aggressive initial conditions %%%%
        x_initial = [1 1 2]';
        v_initial = [2 2 2]';
        roll_initial = pi/4;
        pitch_initial = pi/4;
        yaw_initial = pi/4;
        omega_initial = [1 1 1]';
        f_initial = 0.1;
        f_dot_initial = 0;
     
    otherwise
        error('Unknown IC_CASE value. Please check among the available cases')
end

%% Run ODE
%NO NEED to change em now
rpy_initial = [roll_initial pitch_initial yaw_initial]';
%initialConditions = zeros(14,1); %redundant but kept
initialConditions = [x_initial;v_initial;rpy_initial;omega_initial;f_initial;f_dot_initial];

options = odeset('RelTol',1e-9,'AbsTol',1e-15);
[t,state] = ode45(@(t,state) dfl_approximated_ode(t,state),Tspan,initialConditions,options);

%% Results
x = state(:,1);
y = state(:,2);
z = state(:,3);
yaw = state(:,9);

xd = zeros(length(t),1);
yd = zeros(length(t),1);
zd = zeros(length(t),1);
yawd = zeros(length(t),1);

for i = 1:length(t)
    [xd_i, ~, ~, ~, ~, yawd_i, ~, ~] = reference_planner(t(i));
    xd(i)   = xd_i(1);
    yd(i)   = xd_i(2);
    zd(i)   = xd_i(3);
    yawd(i) = yawd_i;
end

figure(1);plot(t,x,t,xd);legend('x','x_d');xlabel('t [sec]');ylabel('x [m]');title('Position: x(t) and x_d(t)');
figure(2);plot(t,y,t,yd);legend('y','y_d');xlabel('t [sec]');ylabel('y [m]');title('Position: y(t) and y_d(t)');
figure(3);plot(t,z,t,zd);legend('z','z_d');xlabel('t [sec]');ylabel('z [m]');title('Position: z(t) and z_d(t)');
figure(4);plot(t,yaw,t,yawd);legend('yaw','yaw_d');xlabel('t [sec]');ylabel('yaw [m]');title('Position: yaw(t) and yaw_d(t)');
