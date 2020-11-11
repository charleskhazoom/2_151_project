function make_plots(tspan,z_out,u_out,dz_out,ball_alongPlate,accel,p)

    %% plot Energy
    E = energy_arm(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)'); title('Energy Over Time');
    %% plot controls
    figure(2);
    plot(tspan,u_out);
    xlabel('Time (s)'); ylabel('Inputs');
    legend('Force Cart','Torque Joint1', 'Torque Joint2', 'Torque Joint3');

%     %% plot results
    rE = zeros(2,length(tspan));
    vE = zeros(2,length(tspan));
    for i = 1:length(tspan)
        rE(:,i) = position_endEffector(z_out(:,i),p);
        vE(:,i) = velocity_endEffector(z_out(:,i),p);
    end
    
    figure(3); clf;
    plot(tspan,rE(1,:),'r','LineWidth',2)
    hold on
    plot(tspan,rE(2,:),'b','LineWidth',2)
    xlabel('Time (s)'); ylabel('Position (m)'); legend({'x','x_d','y','y_d'});
    title('End Effector position vs desired');

    figure(4); clf;
    plot(tspan,vE(1,:),'r','LineWidth',2)
    hold on
    plot(tspan,vE(2,:),'b','LineWidth',2)
    xlabel('Time (s)'); ylabel('Velocity (m)'); legend({'vel_x','vel_y'});
    title('End Effector velocity vs desired');

    figure(5)
    plot(tspan,z_out(2:4,:)*180/pi)
    legend('q1','q2','q2');
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    title('Joint Angles over trajectory');
    
    figure(6)
    plot(tspan,z_out(6:8,:)*180/pi)
    legend('q1dot','q2dot','q3dot');
    xlabel('Time (s)');
    ylabel('Angular Velocity (deg/sec)');
    title('Joint Velocities over trajectory');
    
    figure(8)
    theta=      180/pi*z_out(2,:)+180/pi*z_out(3,:)-90*ones(1,length(z_out))+180/pi*z_out(4,:);
    thetadot=   180/pi*dz_out(2,:)+180/pi*dz_out(3,:)+180/pi*dz_out(4,:);
    thetadotdot=180/pi*dz_out(6,:)+180/pi*dz_out(7,:)+180/pi*dz_out(8,:); 
    plot(tspan,theta,'b',tspan,thetadot,'g',tspan,thetadotdot,'r')
    legend('angle','rotation','acceleration');
    xlabel('Time (s)');
    ylabel('[deg][deg/sec][deg/sec^2]');
    title('Plate in world frame');
    
    figure(9)
    plot(tspan,ball_alongPlate(1:3,:));
    legend('x','xdot','xdotdot');
    xlabel('Time (s)');
    ylabel('[m][m/s][m/s^2]');
    title('Ball along plate');
    
    figure(10)
    plot(tspan,accel(1:2,:));
    legend('x','y');
    xlabel('Time (s)');
    ylabel('[m/s^2]');
%     axis([0 tf -2 2]);
    title('accleration of plate world frame');
    