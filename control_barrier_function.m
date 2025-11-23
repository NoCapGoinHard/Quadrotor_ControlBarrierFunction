function v=control_barrier_function(t,state,v_star) 

global m g delta delta1 mu2 mu3 mu4

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
%h=(x-obs)'*(x-obs);

% if norm(z)<=delta
%     t
%     error('Collision has happened');
if norm(x-obs)>delta1
    v=v_star;
else
    proj=(z*z')/(z'*z);
    v = obs_4dot -(2*z_dot+mu2*z_2dot+mu3*z_3dot)/mu4 -(mu3*z_dot'*z_2dot+mu4*z_dot'*z_3dot)*(z+mu2*z_dot+mu3*z_2dot+mu4*z_3dot)/(mu4*h) +(eye(3)-proj)*v_star(1:3);

end
