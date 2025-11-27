function v=control_barrier_function(t,state,v_star) 

global next_print_t m g delta delta1 mu1 mu2 mu3 mu4

% quadrotor states
phi = state(7); theta = state(8); psi = state(9);
p = state(10); q = state(11); r = state(12);
f = state(13); f_dot = state(14);
cphi = cos(phi); sphi = sin(phi); ctheta = cos(theta); stheta = sin(theta); cpsi = cos(psi); spsi = sin(psi); ttheta = tan(theta);

% time derivative of rpy
phi_dot = p+sphi*ttheta*q+cphi*ttheta*r;
theta_dot = cphi*q-sphi*r;
psi_dot = sphi/ctheta*q+cphi/ctheta*r;

% position, velocity, accel ecc of the quadrotor
x = state(1:3); 
x_dot = state(4:6);

x_2dot = -f/m*[ cpsi*stheta*cphi+spsi*sphi;
                spsi*stheta*cphi-sphi*cpsi;
                ctheta*cphi] + [0 0 g]';

x_3dot = -f/m*[-sphi*stheta*cpsi*phi_dot+cphi*ctheta*cpsi*theta_dot-cphi*stheta*spsi*psi_dot+...
                cphi*spsi*phi_dot + sphi*cpsi*psi_dot;
                -sphi*stheta*spsi*phi_dot+cphi*ctheta*spsi*theta_dot+cphi*stheta*cpsi*psi_dot+...
                -cphi*cpsi*phi_dot+sphi*spsi*psi_dot;
                -sphi*ctheta*phi_dot-cphi*stheta*theta_dot]+...
         -f_dot/m*[cpsi*stheta*sphi+sphi*spsi;
                   cphi*stheta*spsi-sphi*cpsi;
                   ctheta*cphi];

% position, velocity, accel ecc of the obstacle
[obs,obs_dot,obs_2dot,obs_3dot,obs_4dot]=obstacle(t);
z=x-obs; z_dot=x_dot-obs_dot; z_2dot=x_2dot-obs_2dot; z_3dot=x_3dot-obs_3dot;


% 3th order barrier function
h=z'*(z+mu2*z_dot+mu3*z_2dot+mu4*z_3dot);
h_dot_star=z_dot'*(z+mu2*z_dot+mu3*z_2dot+mu4*z_3dot)+z'*(z_dot+mu2*z_2dot+mu3*z_3dot+mu4*(v_star-obs_4dot));
v_cbf = obs_4dot - ((mu1*z_dot+mu2*z_2dot+mu3*z_3dot)/mu4) - ((mu3*z_dot'*z_2dot+mu4*z_dot'*z_3dot)*(z+mu2*z_dot+mu3*z_2dot+mu4*z_3dot)/(mu4*h)) + (eye(3)-(z*z')/(z'*z))*v_star(1:3);

%if h<=delta
if norm(z) <= delta
    error(['Collision has happened at t = ', num2str(t)]);
end

% the interpolation coefficient "coeff" is introduced in place of 
% the if statement "(h<=delta1)&&(h_dot_star<=0)" in order to 
% smooth the control law, otherwise ode45 crushes because the
% discontinuous controller is too aggressive 

% "a" is equivalent to check "h<=delta1"
%term1 = 2*(1-tanh(1000*(h-delta1)));
%term1=exp(-1000*(h-delta1));
check1=min(exp(-1000*(h-delta1)) , 1 ); % a=1 if h<=delta1, a=0 if h>delta1


% "b" is equivalent to check "h_dot_star<=0"
check2=min((exp(-h_dot_star/0.001)),1); % b=1 if h_dot_star<=0, b=0 if h_dot_star>0

coeff=check1*check2; 

% if coeff=0 then v=v_star is active, if coeff=1 then v=v_cbf 
v = coeff*v_cbf + (1-coeff)*v_star;


if t >= next_print_t
    disp(['t = ', num2str(t), ', h = ', num2str(h), ', h_dot_star = ', num2str(h_dot_star), ', check1 = ', num2str(check1),  ', check2 = ', num2str(check2),  ', check1*check2 = ', num2str(coeff), '.']);
    next_print_t = next_print_t + 0.01;
end

end
