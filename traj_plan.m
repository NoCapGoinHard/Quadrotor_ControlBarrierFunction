function [pd, pd_dot, pd_2dot, pd_3dot, pd_4dot]=traj_plan(t)
% generates a reference trajectory to follow
global traj_type R w
switch traj_type
    case 'line'
        pd      = [t; 0; 0]; pd_dot  = [1; 0; 0]; pd_2dot = [0; 0; 0]; pd_3dot = [0; 0; 0]; pd_4dot = [0; 0; 0];
    case 'circle'
        pd      = [       R*cos(w*t);        R*sin(w*t); 0];
        pd_dot  = [    -R*w*sin(w*t);      R*w*cos(w*t); 0];
        pd_2dot = [-R*(w^2)*cos(w*t); -R*(w^2)*sin(w*t); 0];
        pd_3dot = [ R*(w^3)*sin(w*t); -R*(w^3)*cos(w*t); 0];
        pd_4dot = [ R*(w^4)*cos(w*t);  R*(w^4)*sin(w*t); 0];
    case 'square'
        if t<=10
            pd      = [t; 0; 0]; pd_dot  = [1; 0; 0]; pd_2dot = [0; 0; 0]; pd_3dot = [0; 0; 0]; pd_4dot = [0; 0; 0];
        elseif t<=20
            pd      = [10; 10-t; 0]; pd_dot  = [0; -1; 0]; pd_2dot = [0; 0; 0]; pd_3dot = [0; 0; 0]; pd_4dot = [0; 0; 0];
        elseif t<=30
            pd      = [30-t; -10; 0]; pd_dot  = [-1; 0; 0]; pd_2dot = [0; 0; 0]; pd_3dot = [0; 0; 0]; pd_4dot = [0; 0; 0];
        elseif t<=40
            pd      = [0; -40+t; 0]; pd_dot  = [0; 1; 0]; pd_2dot = [0; 0; 0]; pd_3dot = [0; 0; 0]; pd_4dot = [0; 0; 0];
        end
    otherwise
        error('Please select traj_type among the available values');
end
end