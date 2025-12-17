function animations(t,x,y,xd,yd,xobs,yobs,d_min,selection_vector)
global traj_type M delta1 delta

switch traj_type
    case 'line'
        x_limit = [-1 21];
        y_limit = [-5 5];
    case 'circle'
        x_limit = [-7 7];
        y_limit = [-7 7];
    case 'square'
        x_limit = [-10 10];
        y_limit = [-10 10];
    otherwise
        error('Please select traj_type among the available values');
end


%% plot robot and obstacles trajectories (with time)
if selection_vector(1)~=0
    step_size=1000;
    num_loops = 2;
    tail_length = 100;
    N = length(t);
    indices = 1:step_size:N;

    figure(1); hold on; grid on; axis equal;
    xlabel('x [m]'); ylabel('y [m]');
    title('Robot trajectory with multiple obstacles');
    axis([x_limit(1), x_limit(2), y_limit(1), y_limit(2)]);
    h_ref_line = plot(xd, yd, '--', 'Color', [1 0.4 0.2], 'LineWidth', 1.5);
    h_desired_marker = plot(xd(1), yd(1), 'd', 'MarkerSize', 8, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
    h_robot_path = plot(x(1), y(1), 'Color', [0 0 0.55], 'LineWidth', 1.8);
    h_obs_marker = zeros(M, 1);
    for j = 1:M
        h_obs_marker(j) = plot(xobs(1,j), yobs(1,j), 's', 'MarkerSize', 8, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
    end
    h_robot_marker = plot(x(1), y(1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0 0 0.55], 'MarkerEdgeColor', [0 0 0.55]);
    h_time = text(x(1), y(1) + 0.2, ['t = ' num2str(t(1), '%.2f') ' s'], 'FontSize', 10);
    legend([h_robot_marker, h_robot_path, h_ref_line, h_desired_marker, h_obs_marker(1)], ...
        'Robot Position', 'Robot Path', 'Reference Path', 'Desired Position', 'Obstacles', ...
        'Location', 'best');
    for loop_count = 1:num_loops
        t_prev = t(1);
        set(h_robot_path, 'XData', x(1), 'YData', y(1));
        for idx = indices(2:end)
            start_idx = max(1, idx - tail_length);
            dt = t(idx) - t_prev;
            set(h_robot_path, 'XData', x(start_idx:idx), 'YData', y(start_idx:idx));
            set(h_desired_marker, 'XData', xd(idx), 'YData', yd(idx));
            set(h_robot_marker, 'XData', x(idx), 'YData', y(idx));
            for j = 1:M
                set(h_obs_marker(j), 'XData', xobs(idx, j), 'YData', yobs(idx, j));
            end
            set(h_time, 'String', ['t = ' num2str(t(idx), '%.2f') ' s']);
            set(h_time, 'Position', [x(idx), y(idx) + 0.2]);
            drawnow;
            pause(dt);
            t_prev = t(idx);
        end
        pause(0.5);
    end
    hold off;
end

%% plot robot and obstacles paths (no time)
if selection_vector(2)~=0
    figure(1); hold on; grid on; axis equal;
    plot(x, y, 'Color', [0 0 0.55], 'LineWidth', 1.8);
    plot(xd, yd, '--', 'Color', [1 0.4 0.2], 'LineWidth', 1.5);
    for j = 1:M
        plot(xobs(:,j), yobs(:,j), '--','Color',[0 0 0], 'LineWidth', 2);
    end
    legend('Robot', 'Reference', 'Obstacles');
    xlabel('x [m]'); ylabel('y [m]');
    title('Robot path with multiple obstacles');
    axis([x_limit(1), x_limit(2), y_limit(1), y_limit(2)]);
    %plot(0, -5, 'ko', 'MarkerFaceColor', 'k', 'HandleVisibility', 'off')
    %plot(9, 0, 'ko', 'MarkerFaceColor', 'k', 'HandleVisibility', 'off')
end

%% plot minimum distance
if selection_vector(3)~=0
    figure(2); hold on; grid on;
    plot(t, d_min, 'b', 'LineWidth', 1.8);
    yline(delta,  'r--', 'LineWidth', 2);
    yline(delta1, 'g--', 'LineWidth', 2);
    xlabel('time [s]');
    ylabel('min distance [m]');
    legend('minimum distance', 'collision threshold \delta', 'CBF threshold \delta_1');
    title('Minimum distance from robot to closest obstacle');
end