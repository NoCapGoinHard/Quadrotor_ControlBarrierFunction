function out = dfl_approximated_ode(t,state)

% linear controller to control the feedback linearized model of the quadrotor
% this is the controller when no obstacle in detected
v_star = linear_controller(t,state);

% control barrier function
% for now, it is assumed the obstacle position, velocity ecc is fully known
% in future, an estimation has to be done here 
% notice: yaw (4th entry of v) is not affected by cbf
v=[control_barrier_function(t,state,v_star(1:3));v_star(4)];

% dynamic feedback linearization control
utilde = dynamic_compensator(state,v);
dzeta = [state(14);utilde(1)];

% quadrotor
u = [state(13);utilde(2:4)];
dquadrotor_model = quadrotor_model(state,u);

% ode output
out = [dquadrotor_model;dzeta];

end

