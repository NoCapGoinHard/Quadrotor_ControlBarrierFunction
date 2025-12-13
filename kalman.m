function [obs_dot, obs_2dot, obs_3dot, obs_4dot]=kalman(obs)
global T_s P_obs V obs_est obs_dot_est obs_2dot_est obs_3dot_est obs_4dot_est M

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