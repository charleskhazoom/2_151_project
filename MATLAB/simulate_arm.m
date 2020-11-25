function simulate_arm()
    % simulate_arm: calculates states and inputs over a time span and then
    % plots the energy and kinematics of the system. Also animates the motion
    % of the system
    
    % INPUTS
    % N/A

    % OUTPUTS
    % N/A
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
    
    %k = 1/2; %inertia of cylinder
    k = 2/5; %inertia of ball
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
    
    p_cup_initial = [-0.8,0.5]';
    q0 = eval(invKin_arm(p_cup_initial,p,[0,0]'));
    z0 = [q0;0;0;0;0;0;0]; %8+2[ball: x xdot] last 2 entries for the ball
    
    p_cup_final = [.5,.3]'; % need to specify orientation of last link in the world frame too!
    qf = eval(invKin_arm(p_cup_final,p,q0(2:3)));
    zf = [qf;0;0;0;0;0;0];
    
    z_out = zeros(10,num_step);
    z_out(:,1) = z0;
    ball_alongPlate = zeros(3,num_step);
    ball_alongPlate(:,1) = [0;0;0];
    
    %% choose control law
%     ctrl_law_str = 'joint_space_fb_lin';
%     ctrl_law_str = 'operational_space_fb_lin';
    ctrl_law_str = 'standard_lqr';

    control_law = get_controller(zf(1:8),p_estim,ctrl_law_str); % outputs function handle to be used during Euler integration (for loop below)

% to design a new control law:
% 1) create function in external file which takes as argument
% t,z,u,p,p_ctrl (see example control_law_feedback_linearization.m)
% 2) edit get_controller, where the design parameters are created and
% passed to the anonymous function which is the output of control law.
    %%
    
    for i=1:num_step
        % Compute Controls % HERE WE CAN CHANGE THE CONTROL LAW TO BE ANYTHING

        u_out(:,i) = control_law(tspan(i),z_out(1:8,i));
%         u_out(:,i) = zeros(4,1);
        dz = dynamics(tspan(i), z_out(:,i),u_out(:,i), p);
        %z_out(3,i) = joint_limit_constraint(z_out(:,i),p);
        z_out(:,i+1) = z_out(:,i) + dz*dt;
        dz_out(:,i)=dz;
        
        % Position update
        z_out(1:4,i+1) = z_out(1:4,i) + z_out(5:8,i+1)*dt;
        z_out(9,i+1) = z_out(9,i) + z_out(10,i+1)*dt;
        
        % Ball Simulation
        theta=z_out(2,i)+z_out(3,i)-90/180*pi+z_out(4,i);
        accel(:,i) = acceleration_endEffector([z_out(1:8,i);dz_out(5:8,i)],p);
        a_x_plate = accel(1,i)*cos(-theta) - accel(2,i)*sin(-theta); %rotate into plate frame
        %a_y_plate = accel(1,i)*sin(-theta) + accel(2,i)*cos(-theta); 
        
        a_b_plate = (1/(1+k))*(g*sin(-theta)-a_x_plate); %
        
        ball_alongPlate(3,i+1)=a_b_plate;
        ball_alongPlate(2,i+1)=ball_alongPlate(2,i) + a_b_plate*dt;
        ball_alongPlate(1,i+1)=ball_alongPlate(1,i) + ball_alongPlate(2,i+1)*dt;

    end
    z_out=z_out(:,1:end-1);
    ball_alongPlate=ball_alongPlate(:,1:end-1);
    
%% Look at Results
    make_plots(tspan,z_out,u_out,dz_out,ball_alongPlate,accel,p);

    rE = zeros(2,length(tspan));
    vE = zeros(2,length(tspan));
    for i = 1:length(tspan)
        rE(:,i) = position_endEffector(z_out(:,i),p);
        vE(:,i) = velocity_endEffector(z_out(:,i),p);
    end

    
    %% Animate Solution
    figure(7); clf;
    theta=z_out(2,:)+z_out(3,:)-90/180*pi+z_out(4,:);
    hold on
    animateSol(tspan, z_out,p,ball_alongPlate,rE, theta, p_cup_initial, p_cup_final);
    
end

function animateSol(tspan, x, p, ballX,rEE, theta, start_pos, final_pos)
    % Prepare plot handles
    hold on
    h_ground = plot(linspace(-5,5,100),-3*.0254*ones(1,100),'k','LineWidth',2);
    h_carBase = plot([0],[0],'k','LineWidth',2);
    h_carTop = plot([0],[0],'k','LineWidth',2);
    h_carLSide = plot([0],[0],'k','LineWidth',2);
    h_carRSide = plot([0],[0],'k','LineWidth',2);
    h_link1 = plot([0],[0],'LineWidth',2);
    h_link2 = plot([0],[0],'LineWidth',2);
    h_link3 = plot([0],[0],'LineWidth',2);
    
    start = plot(start_pos(1),start_pos(2),'gx');
    final = plot(final_pos(1),final_pos(2),'rx');
    
    ballPlot= plot([0],[0],'LineWidth',2);
    r=.1; %ball's radius is arbitrary
    
    xlabel('x'); ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    axis([-1 1 -1 1]);

    %Step through and update animation
    for i = 1:length(tspan)
        % skip frame.
        if mod(i,50)
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
        
        %ball=ballX(1,i); %check if it works with manual simulation of ball
        ball=z(9); %using state for ball
        set(ballPlot,'XData',[rEE(1,i)+ball*cosd(-theta(i))  rEE(1,i)+ball*cosd(-theta(i))+r*sind(-theta(i))]);
        set(ballPlot,'YData',[rEE(2,i)-ball*sind(-theta(i))  rEE(2,i)-ball*sind(-theta(i))+r*cosd(-theta(i))]);

        pause(.1)
    end
end