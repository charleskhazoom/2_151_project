function make_plots(tspan, z_out, u_out, dz_out, z_hat_out, dz_hat_out, ball_alongPlate, accel, p, use_obsv, zf)
% make_plots: plot results from simulation - configurations, velocities,
% accelerations, observer behavior, system energy
%
% INPUTS
% tspan: timespan of simulation
% z_out: states over time
% u_out: control input over time
% dz_out: state dynamics over time
% z_hat_out: state estimates over time
% dz_hat_out: esimate dynamics over time
% ball_alongPlate: state of ball on platform
% accel: ball acceleration along platform
% p: parameters
% use_obsv: boolean to determine if observer used to estimate states
% zf: final state

set(0, 'defaultfigurecolor', [1 1 1])

%% energy of arm over time
    E = energy_arm(z_out, p);
    figure(1); clf
    plot(tspan, E);
    xlabel('Time (s)'); ylabel('Energy (J)');
    title('System Energy');
    
%% control input over time
    figure(2);
    subplot(211)
    plot(tspan, u_out(1,:));
    ylabel('Cart Force (N)')
    
    subplot(212)
    plot(tspan, u_out(2:end,:));
    ylabel('Torque (Nm)')
    legend('Torque - Joint1', 'Torque - Joint2', 'Torque - Joint3');

    xlabel('Time (s)'); 
%     title('Control Input')
%     legend('Force - Cart', 'Torque - Joint1', 'Torque - Joint2', 'Torque - Joint3');

%% end effector kinematics
    rE = zeros(2, length(tspan));
    vE = zeros(2, length(tspan));
    
    for i = 1:length(tspan)
        rE(:, i) = position_endEffector(z_out(:, i), p);
        vE(:, i) = velocity_endEffector(z_out(:, i), p);
    end

    % position
    figure(3); clf;
    hold on
    plot(tspan, rE(1, :), 'r', 'LineWidth', 2)
    plot(tspan, rE(2, :), 'b', 'LineWidth', 2)
    xlabel('Time (s)'); ylabel('Position (m)');
    legend({'x', 'y'})
    title('End Effector Position');

    % velocity
    figure(4); clf;
    hold on
    plot(tspan, vE(1, :), 'r', 'LineWidth', 2)
    plot(tspan, vE(2, :), 'b', 'LineWidth', 2)
    xlabel('Time (s)'); ylabel('Velocity (m/s)');
    legend({'vel_x', 'vel_y'});
    title('End Effector Velocity');

%% State variables
    % joint angles
    figure(5)
    subplot(211)
    plot(tspan, z_out(2:4, :)*180/pi)
    legend('q_1', 'q_2', 'q_3');
    ylabel('Joint Angles (deg)');
    
    subplot(212)
    plot(tspan, z_out(1, :))
    xlabel('Time (s)');
    ylabel('Cart Position (m)');
    
    % joint velocities
    figure(6)
    plot(tspan, z_out(6:8, :)*180/pi)
    legend({'$\dot{\theta_1}$', '$\dot{\theta_2}$', '$\dot{\theta_3}$'}, 'Interpreter', 'latex');
    xlabel('Time (s)');
    ylabel('Angular Velocity (deg/sec)');
    title('Joint Velocities');
    
    % plate behavior
    figure(8)
    theta =       180/pi*z_out(2, :) + 180/pi*z_out(3, :) - 90*ones(1, length(z_out)) + 180/pi*z_out(4, :);
    thetadot =    180/pi*dz_out(2, :) + 180/pi*dz_out(3, :) + 180/pi*dz_out(4, :);
    thetadotdot = 180/pi*dz_out(6, :) + 180/pi*dz_out(7, :) + 180/pi*dz_out(8, :); 
    plot(tspan, theta, 'b', tspan, thetadot, 'g', tspan, thetadotdot, 'r')
    legend({'Angle (deg)' , 'Velocity (deg/sec)','Acceleration (deg/sec^2)'});
    xlabel('Time (s)');
    ylabel('Plate Kinematics');
    title('Plate Behavior in World Frame');
    
    % ball behavior
    figure(9)
    plot(tspan, ball_alongPlate(1:3, :));
    legend({'$x$', '$\dot{x}$', '$\ddot{x}$'}, 'Interpreter','latex');
    xlabel('Time (s)');
    ylabel('Ball Kinematics');
    title('Ball along plate');
    
    % plate acceleration
    figure(10)
    plot(tspan, accel(1:2, :));
    legend('x', 'y');
    xlabel('Time (s)');
    ylabel('Acceleration (m/secc^2)');
    title('Accleration of Plate in World Frame');

%% Observer behavior
    str_list = {'$x$ (m)', '$q_1$ (rad)', '$q_2$ (rad)', '$q_3$ (rad)', ...
        '$\dot{x}$ (m/s)', '$\dot{q}_1$ (m/s)', '$\dot{q}_2$ (rad/s)', '$\dot{q}_3$ (rad/s)', ...
        '$x_{b}$ (m)', '$\dot{x}_{b}$ (m/s)'};
    
    figure(11);
    clf;
    for k = 1:10
        subplot(5, 2, k)
        hold on;
        plot(tspan, z_out(k, :)', '-', 'linewidth', 1.5);

        if use_obsv
            plot(tspan, z_hat_out(k, 1:end - 1)', '--', 'linewidth', 1.8);
        end

        plot([tspan(1) tspan(end)], [zf(k) zf(k)], '-.k', 'linewidth', 1.8);
        ylabel(str_list{k}, 'Interpreter', 'latex', 'fontsize', 20)

        if k == 9 || k == 10
            xlabel('Time (s)', 'Interpreter', 'latex', 'fontsize', 20);
        end

        if k == 1
            if use_obsv
                legend({'True State', 'Observer State', 'Desired Final State'}, 'fontsize', 10, 'orientation', 'horizontal')
            else
                legend({'True State', 'Desired Final State'}, 'fontsize', 10, 'orientation', 'horizontal')
            end
        end
    end
end