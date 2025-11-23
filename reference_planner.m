function [xd, xd_dot, xd_2dot, xd_3dot, xd_4dot, yawd, yawd_dot, yawd_2dot] = reference_planner(t)
%REFERENCE_PLANNER  Desired trajectory and derivatives at time t.
%
% Outputs:
%   xd        : desired position [3x1]
%   xd_dot    : first time derivative of xd
%   xd_2dot   : second time derivative of xd
%   xd_3dot   : third time derivative of xd
%   xd_4dot   : fourth time derivative of xd
%   yawd      : desired yaw (scalar)
%   yawd_dot  : first derivative of yawd
%   yawd_2dot : second derivative of yawd

    % Shared configuration (set in SIMULATION_DFL)
    global REF_CASE o

    switch REF_CASE

        case 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1) HEIGHT CONTROL (HOVER) at (0, 0, 3)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            xd      = [t; 0; 0];
            xd_dot  = [1; 0; 0];
            xd_2dot = [0; 0; 0];
            xd_3dot = [0; 0; 0];
            xd_4dot = [0; 0; 0];

            yawd      = 0;
            yawd_dot  = 0;
            yawd_2dot = 0;

        case 2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 2) CIRCUMFERENCE in the XY plane
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Position
            xd      = [cos(o*t);
                       sin(o*t);
                       0];
            % Derivatives
            xd_dot  = [-o*sin(o*t);
                        o*cos(o*t);
                        0];
            xd_2dot = [-o^2*cos(o*t);
                       -o^2*sin(o*t);
                        0];
            xd_3dot = [ o^3*sin(o*t);
                       -o^3*cos(o*t);
                        0];
            xd_4dot = [ o^4*cos(o*t);
                        o^4*sin(o*t);
                        0];

            yawd      = 0;
            yawd_dot  = 0;
            yawd_2dot = 0;

        case 3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 3) ABCD PIECEWISE 9th-ORDER POLYNOMIAL TRAJECTORY
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            A = [0; 0; 0];
            B = [1; 0; 0];
            C = [1; 1; 1];
            D = [0; 1; 1];

            % Duration of each segment
            l = 1.5;  % time to go from A to B (and each other segment)

            % Matrix for 9th-order polynomial boundary conditions
            M = [    l^9        l^8          l^7            l^6          l^5;
                     9*l^8      8*l^7        7*l^6          6*l^5        5*l^4;
                     72*l^7     56*l^6       42*l^5         30*l^4       20*l^3;
                     504*l^6    336*l^5      210*l^4        120*l^3      60*l^2;
                     3024*l^5   1680*l^4     840*l^3        360*l^2      120*l      ];

            % Polynomial coefficients for each segment and each axis
            a_xAB = M \ [B(1)-A(1); 0; 0; 0; 0];
            a_yAB = M \ [B(2)-A(2); 0; 0; 0; 0];
            a_zAB = M \ [B(3)-A(3); 0; 0; 0; 0];

            a_xBC = M \ [C(1)-B(1); 0; 0; 0; 0];
            a_yBC = M \ [C(2)-B(2); 0; 0; 0; 0];
            a_zBC = M \ [C(3)-B(3); 0; 0; 0; 0];

            a_xCD = M \ [D(1)-C(1); 0; 0; 0; 0];
            a_yCD = M \ [D(2)-C(2); 0; 0; 0; 0];
            a_zCD = M \ [D(3)-C(3); 0; 0; 0; 0];

            a_xDA = M \ [A(1)-D(1); 0; 0; 0; 0];
            a_yDA = M \ [A(2)-D(2); 0; 0; 0; 0];
            a_zDA = M \ [A(3)-D(3); 0; 0; 0; 0];

            % Before / after moving around the loop: stay at A
            if t <= 2 || t > 2 + 4*l
                xd      = A;
                xd_dot  = [0; 0; 0];
                xd_2dot = [0; 0; 0];
                xd_3dot = [0; 0; 0];
                xd_4dot = [0; 0; 0];
            else
                % Select segment
                if t > 2 && t <= 2 + l
                    a_x = a_xAB; a_y = a_yAB; a_z = a_zAB;
                    x0  = A;     t0  = 2;
                elseif t > 2 + l && t <= 2 + 2*l
                    a_x = a_xBC; a_y = a_yBC; a_z = a_zBC;
                    x0  = B;     t0  = 2 + l;
                elseif t > 2 + 2*l && t <= 2 + 3*l
                    a_x = a_xCD; a_y = a_yCD; a_z = a_zCD;
                    x0  = C;     t0  = 2 + 2*l;
                elseif t > 2 + 3*l && t <= 2 + 4*l
                    a_x = a_xDA; a_y = a_yDA; a_z = a_zDA;
                    x0  = D;     t0  = 2 + 3*l;
                else
                    % Fallback, should not happen
                    xd      = A;
                    xd_dot  = [0; 0; 0];
                    xd_2dot = [0; 0; 0];
                    xd_3dot = [0; 0; 0];
                    xd_4dot = [0; 0; 0];
                    yawd      = 0;
                    yawd_dot  = 0;
                    yawd_2dot = 0;
                    return;
                end

                tau = t - t0;

                % Position (9th order)
                xd = x0 + [a_x(1)*tau^9 + a_x(2)*tau^8 + a_x(3)*tau^7 + a_x(4)*tau^6 + a_x(5)*tau^5;
                           a_y(1)*tau^9 + a_y(2)*tau^8 + a_y(3)*tau^7 + a_y(4)*tau^6 + a_y(5)*tau^5;
                           a_z(1)*tau^9 + a_z(2)*tau^8 + a_z(3)*tau^7 + a_z(4)*tau^6 + a_z(5)*tau^5];

                % First derivative
                xd_dot = [9*a_x(1)*tau^8 + 8*a_x(2)*tau^7 + 7*a_x(3)*tau^6 + 6*a_x(4)*tau^5 + 5*a_x(5)*tau^4;
                          9*a_y(1)*tau^8 + 8*a_y(2)*tau^7 + 7*a_y(3)*tau^6 + 6*a_y(4)*tau^5 + 5*a_y(5)*tau^4;
                          9*a_z(1)*tau^8 + 8*a_z(2)*tau^7 + 7*a_z(3)*tau^6 + 6*a_z(4)*tau^5 + 5*a_z(5)*tau^4];

                % Second derivative
                xd_2dot = [72*a_x(1)*tau^7 + 56*a_x(2)*tau^6 + 42*a_x(3)*tau^5 + 30*a_x(4)*tau^4 + 20*a_x(5)*tau^3;
                           72*a_y(1)*tau^7 + 56*a_y(2)*tau^6 + 42*a_y(3)*tau^5 + 30*a_y(4)*tau^4 + 20*a_y(5)*tau^3;
                           72*a_z(1)*tau^7 + 56*a_z(2)*tau^6 + 42*a_z(3)*tau^5 + 30*a_z(4)*tau^4 + 20*a_z(5)*tau^3];

                % Third derivative
                xd_3dot = [504*a_x(1)*tau^6 + 336*a_x(2)*tau^5 + 210*a_x(3)*tau^4 + 120*a_x(4)*tau^3 + 60*a_x(5)*tau^2;
                           504*a_y(1)*tau^6 + 336*a_y(2)*tau^5 + 210*a_y(3)*tau^4 + 120*a_y(4)*tau^3 + 60*a_y(5)*tau^2;
                           504*a_z(1)*tau^6 + 336*a_z(2)*tau^5 + 210*a_z(3)*tau^4 + 120*a_z(4)*tau^3 + 60*a_z(5)*tau^2];

                % Fourth derivative
                xd_4dot = [3024*a_x(1)*tau^5 + 1680*a_x(2)*tau^4 + 840*a_x(3)*tau^3 + 360*a_x(4)*tau^2 + 120*a_x(5)*tau;
                           3024*a_y(1)*tau^5 + 1680*a_y(2)*tau^4 + 840*a_y(3)*tau^3 + 360*a_y(4)*tau^2 + 120*a_y(5)*tau;
                           3024*a_z(1)*tau^5 + 1680*a_z(2)*tau^4 + 840*a_z(3)*tau^3 + 360*a_z(4)*tau^2 + 120*a_z(5)*tau];
            end

            yawd      = 0;
            yawd_dot  = 0;
            yawd_2dot = 0;

        otherwise
            error('Unknown REF_CASE value. Please select 1, 2, or 3.');
    end
end
