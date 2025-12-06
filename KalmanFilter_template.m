clc
close all
clear variables
%% simulation settings
R = 1;              % circumpherence radius
w = 1;              % obs. angular velocity
r = 0.15;           % shady trajectory parameter

T = 3*2*pi/w;       % time it takes for a whole circumf. to be completed
n = 1000;           % how many samples in our simulation
T_s = T/n;          % sampling interval (1000 samples in our simulation timespan) 
tspan = 0:T_s:T;    % generated time span vector with sampling times


%% create model of the plant (discrete)

A_obs = [eye(3) T_s.*eye(3) ((T_s^2)/2).*eye(3);
         zeros(3) eye(3) T_s*eye(3);
         zeros(3) zeros(3) eye(3)];

B_obs = zeros([9 1]);
B_noise = zeros([9 1])+1;

C_obs = [eye(3) zeros([3 6])];

D_obs = 0;

% state space model (unspecified sampling time: Ts=-1)
Plant = ss(A_obs, [B_obs B_noise], C_obs, D_obs);

% choose process and output noise covariance
Q_proc = 0.001;
R_meas = [0.01 0 0;
          0 0.01 0;
          0 0 0.01];

%% create Kalman Filter and connecti it to noisy plant

% Kalman Filter and signal name definition
[kalmf,L,P] = kalman(Plant,Q_proc,R_meas);

% specify plant input and output names
%Plant.InputName = {'u', 'w'};                      % scalar valued signals
%Plant.OutputName = {'xt', 'yt', 'zt'};             % vector valued signal (dim. 3)

% create sumblock for injecting measurement noise (output not used here)
%vIn_x = sumblk('x = xt + v(1)');                               
%vIn_y = sumblk('y = yt + v(2)');
%vIn_z = sumblk('z = zt + v(3)');

vIn_x = sumblk('x = x_meas + v(1)');
vIn_y = sumblk('y = y_meas + v(2)');
vIn_z = sumblk('z = z_meas + v(3)');


kalmf.InputName = {'u', 'x', 'y', 'z'};
kalmf.OutputName = {'ye'};


% connect noisy plant to filter (include Plant if needed)
SimModel = connect(vIn_x, vIn_y, vIn_z, kalmf, {'u', 'v(1)', 'v(2)', 'v(3)', 'x_meas', 'y_meas', 'z_meas'}, {'ye'});

%% simulation data
u = zeros([length(tspan) 1]);           % obstacle dynamics not subject to external control
v = zeros([length(tspan) 3]);           % initialize empty measurement noise vector

% generate 'artificial' input to Kalman Filter
y_clean = [R*(-0.9-r*cos(w*tspan)); r*R*sin(w*tspan); zeros([1 n+1])];      % measured trajectory

x_meas = y_clean(1,:)';                                                     % trajectory along x axis
y_meas = y_clean(2,:)';                                                     % trajectory along y axis
z_meas = y_clean(3,:)';                                                     % trajectory along z axis

% white noise generation (based on covariance data)
rng(10,'twister');                                                  % initialize randomness algorithm

w = sqrt(Q_proc)*randn(length(tspan),1);                            % generate process noise vector

v(:,1) = sqrt(R_meas(1,1))*randn(length(tspan),1);                  % generate measurement noise on x
v(:,2) = sqrt(R_meas(2,2))*randn(length(tspan),1);                  % generate measurement noise on y
v(:,3) = sqrt(R_meas(3,3))*randn(length(tspan),1);                  % generate measurement noise on z

out = lsim(SimModel, [u v(:,1) v(:,2) v(:,3) x_meas y_meas z_meas], tspan);

%% compute state prediction x[n+1|n]

state_est = out(:,4:12);                 % select state estimation from simulation output
noisy_y = y_clean'+v;                    % column vector of noisy measured output

% time instants for prediction attempts
n_pred = 10;                            % choose how many prediction instants in simulation
pred_instants = T/n_pred;               % increments in sim. timespan for prediction
pred_tspan = 0:pred_instants:T;         % generate prediction instants in the timespan

custom_state_est = zeros([9 n_pred]);    % generate empty vector for estimations at specific instants
    
% compute state increments using Kalman Gain
state_dot = zeros([9 n_pred]);          % generate empty verctor to fill with increments

for i=1:n_pred                          % fill variation and estimates vectors
    t_pred = i*(n/n_pred)+1;
    custom_state_est(:,i) = state_est(t_pred,:);
    state_dot(:,i) = A_obs*state_est(t_pred,:)' + L*(noisy_y(t_pred,:)' - C_obs*state_est(t_pred,:)');
end

% compute state prediction x[n+1|n]
state_pred = custom_state_est + state_dot .* T_s;

%% predict trajectory at given time instants

n_output = 3;                                       % specify number of measurements (model outputs)
horizon = 5;                                       % set prediction horizon
y_pred = zeros([n_output*horizon n_pred]);          % generate empty vector for trajectory predictions
 
for i=1:n_pred
    
    pred = traj_pred(state_pred(:,i), A_obs, C_obs, horizon, n_output);   % call function to get trajectory prediction

    for j=1:horizon
        y_pred((j-1)*n_output+1:j*n_output,i) = (pred(j,:))';        % put the output of the function in a single column in y_pred
    end
end

%% plot predictions

xt = x_meas;             % true x response
xe = out(:,1);           % filtered x response
x = xt + v(:,1);         % measured x response

x_pred = zeros([horizon n_pred]);       % generate empty vector for x coord. prediction

for i=1:n_pred                          % fill empty vector with predictions (x coord.)
    for j=1:horizon
        x_pred(j,i) = y_pred((j-1)*n_output+1,i);
    end
end

x_pred_tspan = tspan .* 0;              % generate empty tspan vector to plot predictions

for i=1:n_pred-1                        % fill empty vector with predictions (x coord.)
    t_pred = i*(n/n_pred)+1;
    for j=1:horizon
        x_pred_tspan(t_pred:t_pred+horizon-1) = x_pred(:,i);
    end
end


%% compare responses
%{
xt = x_meas;             % true x response
xe = out(:,1);           % filtered x response
x = xt + v(:,1);         % measured x response

yt = y_meas;             % true y response
ye = out(:,2);           % filtered y response              
y = yt + v(:,2);         % measured y response

clf

subplot(311)
plot(tspan,xt,'b',tspan, xe,'r--') 
xlabel('Number of Samples')
ylabel('Output (x)')
title('Kalman Filter Response')
legend('True','Filtered')

subplot(312)
plot(tspan,yt,'b',tspan, ye,'r--') 
xlabel('Number of Samples')
ylabel('Output (y)')
legend('True','Filtered')

subplot(313)
plot(tspan, yt-ye,'g',tspan,xt-xe,'r--')
xlabel('Number of Samples')
ylabel('Error')
legend('y_{true}-y_{filter}','x_{true} - x_{filter}')
%}


%% FUNCTION: predict trajectory given a state prediction x[n+1|n]
function traj_pred = traj_pred(state_pred_col, A, C, N, n_output,Ts)
    
    traj_pred = zeros([N n_output]);

    for i=1:N
        traj_pred(i,:) = (Ts^2)*(A^(i-1)*state_pred_col) ;        % put predicted output at step k+i in a row vector
    end
end