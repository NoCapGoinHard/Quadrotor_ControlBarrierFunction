function obs_all=obs_traj_multi(t)
% generates the obstacle trajectory
global traj_type w R M
obs_all     = zeros(3,M);
switch traj_type
    case 'line'
        % OBSTACLE 1: line
        vel=1;
        alpha=0;
        %obs_all(:,1)=[-vel*t+5*(1+vel); vel*(t-5)*tan(alpha); 0];
        obs_all(:,1)=[5; vel*(-t+5); 0];
        % OBSTACLE 2: static
        obs_all(:,2)     = [9; 0; 0];
        % OBSTACLE 3: parabola
        obs_all(:,3)=[t; t^2-27*t+180; 0];
    case 'circle'
        % OBSTACLE 1: static
        obs_all(:,1)     = [0; -1; 0];
        % OBSTACLE 2: vertical line
        obs_all(:,2)     = [-1; pi - t; 0];
        % OBSTACLE 3: circle
        r = 0.3;
        obs_all(:,3)     = [R*(0 - r*cos(w*t)); 1 + r*R*sin(w*t); 0];
    case 'square'
        % OBSTACLE 1: circle
        r = 3;
        obs_all(:,1)     = [r*cos(-w*t+3.1); r*sin(-w*t+3.1); 0];
        % OBSTACLE 2: parabola
        obs_all(:,2)     = [t^2-30*t+231; 10-t; 0];
        % OBSTACLE 3: line
        obs_all(:,3)     = [-1+5*cos(t); -6-5*cos(t); 0];
    otherwise
        error('Please select traj_type among the available values');
end
end