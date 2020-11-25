function make_plots(tspan, z_out, u_out, dz_out, ball_alongPlate, accel, p)
% make_plots: plot results from simulation - configurations, velocities,
% accelerations
%
% INPUTS
% tspan: timespan of simulation
% z_out: states over time
% dz_out: rate of change of states over time
% ball_alongPlate: state of ball on platform
% accel: ball acceleration along platform
% p: parameters

%% plot energy of arm over time
    E = energy_arm(z_out, p);
    figure(1); clf
    plot(tspan, E);
    xlabel('Time (s)'); ylabel('Energy (J)');
    title('System Energy');
    
%% plot control input over time
    figure(2);
    plot(tspan, u_out);
    xlabel('Time (s)'); ylabel('Inputs');
    title('Control Input')
    legend('Force - Cart', 'Torque - Joint1', 'Torque - Joint2', 'Torque - Joint3');

%% plot results over time
    rE = zeros(2, length(tspan));
    vE = zeros(2, length(tspan));
    
    for i = 1:length(tspan)
        rE(:, i) = position_endEffector(z_out(:, i), p);
        vE(:, i) = velocity_endEffector(z_out(:, i), p);
    end
    
    % end effector position
    figure(3); clf;
    hold on
    plot(tspan, rE(1, :), 'r', 'LineWidth', 2)
    plot(tspan, rE(2, :), 'b', 'LineWidth', 2)
    xlabel('Time (s)'); ylabel('Position (m)');
%     legend({'x', 'x_d_e_s', 'y', 'y_d_e_s'});
    legend({'x', 'y'})
    title('End Effector Position');

    % end effector velocity
    figure(4); clf;
    hold on
    plot(tspan, vE(1, :), 'r', 'LineWidth', 2)
    plot(tspan, vE(2, :), 'b', 'LineWidth', 2)
    xlabel('Time (s)'); ylabel('Velocity (m/s)');
    legend({'vel_x', 'vel_y'});
    title('End Effector Velocity');

    % joint angles
    figure(5)
    plot(tspan, z_out(2:4, :)*180/pi)
    legend('\theta_1', '\theta_2', '\theta_3');
    xlabel('Time (s)'); ylabel('Angle (deg)');
    title('Joint Angles');
    
    % joint velocities
    figure(6)
    plot(tspan, z_out(6:8, :)*180/pi)
    legend({'$\dot{\theta_1}$', '$\dot{\theta_2}$', '$\dot{\theta_3}$'}, 'Interpreter','latex');
    xlabel('Time (s)');
    ylabel('Angular Velocity (deg/sec)');
    title('Joint Velocities');
    
    % plate behavior
    figure(8)
    theta =       180/pi*z_out(2, :) + 180/pi*z_out(3, :) - 90*ones(1, length(z_out)) + 180/pi*z_out(4, :);
    thetadot =    180/pi*dz_out(2, :) + 180/pi*dz_out(3, :) + 180/pi*dz_out(4, :);
    thetadotdot = 180/pi*dz_out(6, :) + 180/pi*dz_out(7, :) + 180/pi*dz_out(8, :); 
    plot(tspan, theta, 'b', tspan, thetadot, 'g', tspan, thetadotdot, 'r')
    legend({'Angle (deg)' , 'Velocity (deg/sec)','Acceleration (deg/sec^2'});
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
end