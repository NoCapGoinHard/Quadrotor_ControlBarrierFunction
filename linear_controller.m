%%%%%%%Tesi
%%%%%

function v = linear_controller(t,state)

%% initializing variables

x = state(1:3); x_dot = state(4:6);
phi = state(7); theta = state(8); psi = state(9);
p = state(10); q = state(11); r = state(12);
omega = state(10:12);
f = state(13); f_dot = state(14);

global m g J
cphi = cos(phi); sphi = sin(phi); ctheta = cos(theta); stheta = sin(theta); cpsi = cos(psi); spsi = sin(psi); ttheta = tan(theta);
Jinv = inv(J);
J1X = Jinv(1,:); J2X = Jinv(2,:); J3X = Jinv(3,:);
cross_omega_Jomega = cross(omega,J*omega);

global c0 c1 c2 c3 o

%% planning

[xd, xd_dot, xd_2dot, xd_3dot, xd_4dot, yawd, yawd_dot, yawd_2dot] = reference_planner(t);


%% time derivative of rpy
phi_dot = p+sphi*ttheta*q+cphi*ttheta*r;
theta_dot = cphi*q-sphi*r;
psi_dot = sphi/ctheta*q+cphi/ctheta*r;

%% computing derivatives

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

yaw = psi;
yaw_dot = psi_dot;

%% errors

e0 = xd-x;
e1 = xd_dot-x_dot;
e2 = xd_2dot-x_2dot;
e3 = xd_3dot-x_3dot;

e1yaw   = yawd_dot-yaw_dot;
e0yaw   = yawd-yaw;

%% output control law

vx = xd_4dot + c3*e3 + c2*e2 + c1*e1 + c0*e0;
vyaw = yawd_2dot + c1*e1yaw + c0*e0yaw;

v = [vx;vyaw];

end
