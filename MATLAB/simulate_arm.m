function simulate_arm()
    clear; clc
    %% Define fixed paramters
    m_cart =50 ;         
    m1 =5; 
    m2 =5;            
    m3 =5;
    h_cart = 6*0.0254;      
    l_cart = 12*0.0254;
    l_1 = 12*0.0254;
    l_2 = 12*0.0254;
    l_3 = 12*0.0254;
    g = 9.81;    
    
    m_ball = .05;
    u=.25;
    %% Parameter vector (real system)
    p   = [m_cart m1 m2 m3 h_cart l_cart l_1 l_2 l_3 g]';        % parameters
    %% Parameter vector (estimated system)

    p_estim = p; % e
    p_estim(1:4)=p_estim(1:4)*1; % assume all masses are overestimated by some percentage
    
   
    %% Setup Dynamic simulation
    dt = 0.001;
    tf = 10; %May have to change if 10 second not enough to complete task
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 
    
    p_cup_initial = [-0.8,0.85]';
    q0 = eval(invKin_arm(p_cup_initial,p,[0,0,0,0]'));
    z0 = [q0;0;0;0;0];
    
    p_cup_final = [1,.5]'; % need to specify orientation of last link in the world frame too!
    qf = eval(invKin_arm(p_cup_final,p,q0));
    zf = [qf;0;0;0;0];
    
    z_out = zeros(8,num_step);
    z_out(:,1) = z0;
    ball_alongPlate = zeros(3,num_step);
    ball_alongPlate(:,1) = [0;0;0];
    
    %% tangent linearization of dynamics about final desired position
    u_equi = Grav_arm(zf,p); % inputs at equilibrium = gravity
    [A_lin, B_lin] =  linearize_dynamics(zf,u_equi,p);
    
    %% choose control law
    ctrl_law_str = 'joint_space_fb_lin';
%     ctrl_law_str = 'operational_space_fb_lin';

    control_law = get_controller(zf,p_estim,ctrl_law_str); % outputs function handle to be used during Euler integration (for loop below)

% to design a new control law:
% 1) create function in external file which takes as argument
% t,z,u,p,p_ctrl (see example control_law_feedback_linearization.m)
% 2) edit get_controller, where the design parameters are created and
% passed to the anonymous function which is the output of control law.
    %%
    
    for i=1:num_step
        % Compute Controls % HERE WE CAN CHANGE THE CONTROL LAW TO BE ANYTHING
        p_ctrl_fl_op_sp =1;% dummy assignation for now;
%         u_out(:,i) = control_law_fb_lin_op_sp(tspan(i), z_out(:,i), p,p_ctrl_fl_op_sp);
        
        u_out(:,i) = control_law(tspan(i),z_out(:,i));
        dz = dynamics(tspan(i), z_out(:,i),u_out(:,i), p);
        %z_out(3,i) = joint_limit_constraint(z_out(:,i),p);
        z_out(:,i+1) = z_out(:,i) + dz*dt;
        dz_out(:,i)=dz;
        
        % Position update
        z_out(1:4,i+1) = z_out(1:4,i) + z_out(5:8,i+1)*dt;
        
        % Ball Simulation
        theta=z_out(2,i)+z_out(3,i)-90/180*pi+z_out(4,i);
        ddtheta =dz_out(6,i)+dz_out(7,i)+dz_out(8,i);
        accel(:,i) = acceleration_endEffector([z_out(:,i);dz_out(5:8,i)],p);
        a_x_plate = accel(1,i)*cos(-theta) - accel(2,i)*sin(-theta);
        a_y_plate = accel(1,i)*sin(-theta) + accel(2,i)*cos(-theta);
        
        Normal = m_ball*g*cos(-theta) - ddtheta*m_ball*ball_alongPlate(1,i)^3 + m_ball*a_y_plate;
        if(abs(u*Normal)>abs(m_ball*g*sin(-theta)))
            friction=m_ball*g*sin(-theta);
        else
            friction=u*Normal;
        end
        a_b_plate = m_ball*g*sin(-theta) - friction - m_ball*a_x_plate;
        
        ball_alongPlate(3,i+1)=a_b_plate;
        ball_alongPlate(2,i+1)=ball_alongPlate(2,i) + a_b_plate*dt;
        ball_alongPlate(1,i+1)=ball_alongPlate(1,i) + ball_alongPlate(2,i+1)*dt;

    end
    z_out=z_out(:,1:end-1);
    ball_alongPlate=ball_alongPlate(:,1:end-1);
    
    %% plot Energy
    E = energy_arm(z_out,p);
    figure(1); clf
    %plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)'); title('Energy Over Time');
    %% plot controls
    figure(2);
    plot(tspan,u_out);
    xlabel('Time (s)'); ylabel('Inputs');
    legend('Force Cart','Torque Joint1', 'Torque Joint2', 'Torque Joint3');

    %% plot results
    rE = zeros(2,length(tspan));
    vE = zeros(2,length(tspan));
    for i = 1:length(tspan)
        rE(:,i) = position_endEffector(z_out(:,i),p);
        vE(:,i) = velocity_endEffector(z_out(:,i),p);
    end
    
    figure(3); clf;
    plot(tspan,rE(1,:),'r','LineWidth',2)
    hold on
%     plot(tspan,p_traj.x_0 + p_traj.r * cos(p_traj.omega*tspan) ,'r--');
    plot(tspan,rE(2,:),'b','LineWidth',2)
%     plot(tspan,p_traj.y_0 + p_traj.r * sin(p_traj.omega*tspan) ,'b--');
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
    axis([0 tf -2 2]);
    title('accleration of plate world frame');
    
    %% Animate Solution
    figure(7); clf;
    hold on
   
    animateSol(tspan, z_out,p,ball_alongPlate,rE, theta);
end



function dz = dynamics(t,z,u,p)
    % Get mass matrix
    A = A_arm(z,p);
     
    % Get b = Q - V(q,qd) - G(q)
    b = b_arm(z,u,p);
    
    % Solve for qdd.
    qdd = A\(b);
    dz = 0*z;
    
    % Form dz
    dz(1:4) = z(5:8);
    dz(5:8) = qdd;
end


function animateSol(tspan, x, p, ballX,rEE, theta)
    % Prepare plot handles
    hold on
    h_ground = plot(linspace(-5,5,100),-0.15*ones(1,100),'k','LineWidth',2);
    h_carBase = plot([0],[0],'k','LineWidth',2);
    h_carTop = plot([0],[0],'k','LineWidth',2);
    h_carLSide = plot([0],[0],'k','LineWidth',2);
    h_carRSide = plot([0],[0],'k','LineWidth',2);
    h_link1 = plot([0],[0],'LineWidth',2);
    h_link2 = plot([0],[0],'LineWidth',2);
    h_link3 = plot([0],[0],'LineWidth',2);
    
    ball= plot([0],[0],'LineWidth',2);
    r=.1; %ball's radius is arbitrary
    
    xlabel('x'); ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    axis([-1 1 -1 1]);

    %Step through and update animation
    for i = 1:length(tspan)
        % skip frame.
        if mod(i,10)
            continue;
        end
        t = tspan(i);
        z = x(:,i); 
        keypoints = keypoints_arm(z,p);

        rA = keypoints(:,1); % center of mass (which is at 0) of cart
        rB = keypoints(:,2); % where link 1 and 2 meet
        rC = keypoints(:,3); % where link 2 and 3 meet
        rD = keypoints(:,4); % where right side of plate(3)
        rE = keypoints(:,5); % where Left side of plate(3)

        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        set(h_carLSide,'XData',[rA(1)-p(6)/2 rA(1)-p(6)/2]);
        set(h_carLSide,'YData',[rA(2)-p(5)/2 rA(2)+p(5)/2]);
        
        set(h_carRSide,'XData',[rA(1)+p(6)/2 rA(1)+p(6)/2]);
        set(h_carRSide,'YData',[rA(2)-p(5)/2 rA(2)+p(5)/2]);
        
        set(h_carBase,'XData',[rA(1)-p(6)/2 rA(1)+p(6)/2]);
        set(h_carBase,'YData',[rA(2)-p(5)/2 rA(2)-p(5)/2]);
        
        set(h_carTop,'XData',[rA(1)-p(6)/2 rA(1)+p(6)/2]);
        set(h_carTop,'YData',[rA(2)+p(5)/2 rA(2)+p(5)/2]);
        
        set(h_link1,'XData',[rA(1) rB(1)]);
        set(h_link1,'YData',[rA(2) rB(2)]);
        
        set(h_link2,'XData',[rB(1) rC(1)]);
        set(h_link2,'YData',[rB(2) rC(2)]);
        
        set(h_link3,'XData',[rE(1) rD(1)]);
        set(h_link3,'YData',[rE(2) rD(2)]);
        
        set(ball,'XData',[rEE(1,i)+ballX(1,i)*cosd(-theta(i))  rEE(1,i)+ballX(1,i)*cosd(-theta(i))+r*sind(-theta(i))]);
        set(ball,'YData',[rEE(2,i)-ballX(1,i)*sind(-theta(i))  rEE(2,i)-ballX(1,i)*sind(-theta(i))+r*cosd(-theta(i))]);

        pause(.01)
    end
end


function [A_lin, B_lin] =  linearize_dynamics(z_equi,u_equi,p)
    n_states = length(z_equi);
    n_inputs = length(u_equi);
    t =0; % dummy assignation of t variable
    dz_equi = dynamics(t,z_equi,u_equi,p);
%     zout0 = zout(:,end);
    
    %statef0
    ep = 1e-6;
    delStateMat = eye(n_states);
    A_lin = zeros(n_states,n_states);
    B_lin = zeros(n_states,n_inputs);

    for i=1:n_states
        z_dev = z_equi + ep*delStateMat(:,i);
        dz_dev = dynamics(t,z_dev,u_equi,p);
        A_lin(:,i) = (dz_dev - dz_equi)/ep;
    end


    del_u_mat= eye(n_inputs);
    for i=1:n_inputs
        u_dev = u_equi+ep*del_u_mat(:,i);
        dz_dev = dynamics(t,z_equi,u_dev,p);
        B_lin(:,i) = (dz_dev - dz_equi)/ep;
    end
end

