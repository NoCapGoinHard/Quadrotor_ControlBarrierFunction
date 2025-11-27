%% Dynamic Feedback Linearization of the quadrotor complete model
%% and Control Barrier Function for obstacle avoidance

%% November 2025
%% Carlo Rugiero, Francesco Maria Germano, Giovanni Pio Cuoco, Matteo Trusiani 
% based on software developed by Paolo Forni for his bachelor degree thesis (June 2012) 
% supervisor: Marilena Vendittelli

clear variables; close all; clc

%% Parameters
global next_print_t Ix Iy Iz g m J c0 c1 c2 c3 REF_CASE delta delta1 mu1 mu2 mu3 mu4

next_print_t = 0.01;

T = 5.5; % total simulation time

% quadrotor parameters
g = 9.81;
% m = 4.34; % Ix = 0.0820; % Iy = 0.0845; % Iz = 0.1377;
m = 4;
Ix = 0.1;
Iy = 0.1;
Iz = 0.1;
J = [Ix 0 0; 0 Iy 0; 0 0 Iz];

% linear controller parameters: Hurwitz polynomial coefficients
% to assign the characteristic polynomial (s-eig)^4
eig=10;
c3 = 4*eig; c2 = 6*eig^2; c1 = 4*eig^3; c0 = 1*eig^4;    %%c3 = 26;c2 = 253;c1 = 1092; c0 = 1764; % to assigns eigenvalues: -6,-6,-7,-7

% control barrier function parameters
delta = 0.05;
delta1 = 0.5;
mu=0.05; % 0.05
q=8; % 8
mu1 = 2;
mu2 = mu*q;
mu3 = 0;
mu4 = mu;

% reference trajectory:
% 1 = hover, 2 = circumference, 3 piecewise polynomial
REF_CASE = 1;

% initial conditions:
initialConditions=zeros(1,14); % [x y z vx vy vz phi theta psi wphi wtheta wpsi f f_dot]
initialConditions(13)=0.1; % f cannot be zero otherwise the decoupling matrix will become singular

%% Run ODE

Tspan = [0 T]; 
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
