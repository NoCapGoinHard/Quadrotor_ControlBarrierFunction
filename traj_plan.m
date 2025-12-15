function [pd, pd_dot, pd_2dot, pd_3dot, pd_4dot]=traj_plan(t)
% generates a reference trajectory to follow
global traj_type
switch traj_type
    case 'line'
        pd      = [t; 0; 0]; pd_dot  = [1; 0; 0]; pd_2dot = [0; 0; 0]; pd_3dot = [0; 0; 0]; pd_4dot = [0; 0; 0];
    case 'circle'
        R = 5; % radious of circle
        pd      = [ R*cos(t);  R*sin(t); 0];
        pd_dot  = [-R*sin(t);  R*cos(t); 0];
        pd_2dot = [-R*cos(t); -R*sin(t); 0];
        pd_3dot = [ R*sin(t); -R*cos(t); 0];
        pd_4dot = [ R*cos(t);  R*sin(t); 0];
    case 'square'
        if t<=10
            pd      = [t-5; 5; 0]; pd_dot  = [1; 0; 0]; pd_2dot = [0; 0; 0]; pd_3dot = [0; 0; 0]; pd_4dot = [0; 0; 0];
        elseif t<=20
            pd      = [5; 15-t; 0]; pd_dot  = [0; -1; 0]; pd_2dot = [0; 0; 0]; pd_3dot = [0; 0; 0]; pd_4dot = [0; 0; 0];
        elseif t<=30
            pd      = [25-t; -5; 0]; pd_dot  = [-1; 0; 0]; pd_2dot = [0; 0; 0]; pd_3dot = [0; 0; 0]; pd_4dot = [0; 0; 0];
        elseif t<=40
            pd      = [-5; -35+t; 0]; pd_dot  = [0; 1; 0]; pd_2dot = [0; 0; 0]; pd_3dot = [0; 0; 0]; pd_4dot = [0; 0; 0];
        end
    otherwise
        error('Please select traj_type among the available values');
end
end