function simulate_arm()
% simulate_arm: calculates states and inputs over a time span and then
% plots the energy and kinematics of the system. Also animates the motion
% of the system
%
% INPUTS
% N/A
%
% OUTPUTS
% N/A

    clear; clc;
    
%% Define fixed paramters
    m_cart = 50; % cart mass, kg         
    m1 = 5; % link 1 mass, kg
    m2 = 5; % link 2 mass, kg
    m3 = 5; % link 3 mass, kg
    h_cart = 6*0.0254; % cart height, m
    l_cart = 12*0.0254; % cart length, m
    l_1 = 12*0.0254; % link 1 length, m
    l_2 = 12*0.0254; % link 2 length, m
    l_3 = 12*0.0254; % link 3 length, m
    g = 9.81; % gravitational acceleration, m/s^2
    
    m_ball = .1; % ball mass, kg
    % seems like making mass larger will result in ball moving more
    mu =.25; % coefficient of friction between ball and platform
    
%% Parameter vector (real system)
    p = [m_cart m1 m2 m3 h_cart l_cart l_1 l_2 l_3 g]'; % parameters vector
   
%% Parameter vector (estimated system)
    p_estim = p;
    
    mass_error = 1; % assume all masses are estimated incorrectly by some percentage
    p_estim(1:4) = p_estim(1:4)*mass_error;
   
%% Set up dynamic simulation
    % State vector:
    % [x_cart, th_1, th_2, th_3, vx_cart, omega_1, omega_2, omega_3]'
    
    dt = 0.001; % timestep, sec
    tf = 10; % % final time, sec (change if 10 seconds not enough to complete task)
    
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 
    
    p_cup_initial = [-0.8, 0.5]'; % initial position of ball (x, y), m
    q0 = eval(invKin_arm(p_cup_initial, p , [0, 0]')) % initial configuration
    z0 = [q0; 0; 0; 0; 0]; % initial state
    
    p_cup_final = [0.5, 0.3]'; % need to specify orientation of last link in the world frame too!
    qf = eval(invKin_arm(p_cup_final,p,q0(2:3))) % final configuration
    zf = [qf; 0; 0; 0; 0]; % final state
    
    z_out = zeros(8, num_step); % state evolution over time
    z_out(:, 1) = z0;
    
    dz_out = zeros(8, num_step); % store rate of change of states at each time step
    
    u_out = zeros(4, num_step); % store control input at each time step 
    
    ball_alongPlate = zeros(3, num_step); % ball position over time
    ball_alongPlate(:, 1) = [0; 0; 0];
    accel = zeros(2, num_step); % acceleration of platform
    
%% Choose control law
    ctrl_law_str = 'joint_space_fb_lin';
%     ctrl_law_str = 'operational_space_fb_lin';
%     ctrl_law_str = 'standard_lqr';

    control_law = get_controller(zf, p_estim, ctrl_law_str);
    % 'control_law' is a function handle to be used during Euler integration
    % Function handle will take time and current state and return next control input 
    
    % to design a new control law:
    % 1) create function in external file which takes as argument
    %    t, z, u, p, p_ctrl (see example control_law_feedback_linearization.m)
    % 2) edit get_controller, where the design parameters are created and
    %    passed to the anonymous function which is the output of control law.
    
%% Perform Euler integration
    for i = 1:num_step
        % Compute Controls
        
        % get control input for this timestep
        u_out(:, i) = control_law(tspan(i), z_out(:, i));
        % u_out(:, i) = zeros(4, 1);
    
        % calculate dz, change in state variables
        dz = dynamics(tspan(i), z_out(:, i), u_out(:, i), p);
        % z_out(3, i) = joint_limit_constraint(z_out(:, i), p);
        
        % Forward Euler to calculate next state
        z_out(:, i + 1) = z_out(:, i) + dz*dt;
        dz_out(:, i) = dz;
        
        % Position update
        % new_pos = curr_pos + curr_vel*timestep
        z_out(1:4, i + 1) = z_out(1:4, i) + z_out(5:8, i + 1)*dt;
        
        % Ball Simulation
        theta = z_out(2, i) + z_out(3, i) - 90/180*pi + z_out(4, i); % orientation of platform
        ddtheta = dz_out(6, i) + dz_out(7, i) + dz_out(8, i); % angular velocity of platform
        
        accel(:, i) = acceleration_endEffector([z_out(:, i); dz_out(5:8, i)], p); % acceleration of platform
        a_x_plate = accel(1, i)*cos(-theta) - accel(2, i)*sin(-theta);
        a_y_plate = accel(1, i)*sin(-theta) + accel(2, i)*cos(-theta);
        
        Normal = m_ball*g*cos(-theta) - ddtheta*m_ball*ball_alongPlate(1, i)^3 + m_ball*a_y_plate; % normal force
        % frictional force
        if(abs(mu*Normal) > abs(m_ball*g*sin(-theta)))
            friction = m_ball*g*sin(-theta);
        else
            friction = mu*Normal;
        end
        
        a_b_plate = m_ball*g*sin(-theta) - friction - m_ball*a_x_plate;
        
        % update ball position on platform
        ball_alongPlate(3, i + 1) = a_b_plate;
        ball_alongPlate(2, i + 1) = ball_alongPlate(2, i) + a_b_plate*dt;
        ball_alongPlate(1, i + 1) = ball_alongPlate(1, i) + ball_alongPlate(2, i + 1)*dt;
    end
    
    z_out = z_out(:, 1:end - 1); % final state
    ball_alongPlate = ball_alongPlate(:, 1:end - 1); % final ball position
    
%% Look at Results
    make_plots(tspan, z_out, u_out, dz_out, ball_alongPlate, accel, p);

    rE = zeros(2, length(tspan)); % end effector position
    vE = zeros(2, length(tspan)); % end effector velocity
    
    for i = 1:length(tspan)
        rE(:, i) = position_endEffector(z_out(:, i), p);
        vE(:, i) = velocity_endEffector(z_out(:, i), p);
    end
    
%% Animate Solution
    figure(7); clf;
    theta = z_out(2, :) + z_out(3, :) - 90/180*pi + z_out(4, :); % platform angle
    hold on
    animateSol(tspan, z_out, p, ball_alongPlate, rE, theta, p_cup_initial, p_cup_final);
end

function animateSol(tspan, x, p, ballX, rEE, theta, start_pos, final_pos)
% animateSol: animate robot and ball positions over time using evolution of
% the states and control found above
%
% INPUTS
% tspan: timespan over simulation was run
% x: states
% p: system parameters
% ballX: ball state
% rEE: position of end effector
% theta: platform angle
% start_pos: initial position
% final_pos: goal position
   
    % Prepare plot handles
    hold on
    h_ground = plot(linspace(-5, 5, 100), -3*.0254*ones(1, 100), 'k', 'LineWidth', 2); % ground
    h_carBase = plot([0], [0], 'k', 'LineWidth', 2); % bottom of cart
    h_carTop = plot([0], [0], 'k', 'LineWidth', 2); % top of cart
    h_carLSide = plot([0], [0], 'k', 'LineWidth', 2); % left side of cart
    h_carRSide = plot([0], [0], 'k', 'LineWidth', 2); % right side of cart
    h_link1 = plot([0], [0], 'LineWidth', 2); % link 1
    h_link2 = plot([0], [0], 'LineWidth', 2); % link 2
    h_link3 = plot([0], [0], 'LineWidth', 2); % link 3
    
    start = plot(start_pos(1), start_pos(2), 'gx'); % intial position
    final = plot(final_pos(1), final_pos(2), 'rx'); % final/goal position
    
    ball = plot([0], [0], 'LineWidth', 2); % initial ball position
    r = 0.1; % ball radius, m (arbitrary)
    
    xlabel('x'); ylabel('y');
    h_title = title('t = 0.0 s');
    
    axis equal
    axis([-1 1 -1 1]);

    % Step through and update animation
    for i = 1:length(tspan)
        
        % skip each 50th frame
        if mod(i, 50)
            continue;
        end
        
        t = tspan(i); % time
        z = x(:, i); % state
        keypoints = keypoints_arm(z, p); % defining points  of arm

        rA = keypoints(:, 1); % center of mass (which is at 0) of cart
        rB = keypoints(:, 2); % where link 1 and 2 meet
        rC = keypoints(:, 3); % where link 2 and 3 meet
        rD = keypoints(:, 4); % right side of plate(3)
        rE = keypoints(:, 5); % left side of plate(3)

        set(h_title, 'String', sprintf('t = %.2f', t)); % update title with time
        
        % plot left side of cart
        set(h_carLSide, 'XData', [rA(1) - p(6)/2, rA(1) - p(6)/2]);
        set(h_carLSide, 'YData', [rA(2) - p(5)/2, rA(2) + p(5)/2]);
        
        % plot right side of cart
        set(h_carRSide, 'XData', [rA(1) + p(6)/2, rA(1) + p(6)/2]);
        set(h_carRSide, 'YData', [rA(2) - p(5)/2, rA(2) + p(5)/2]);
        
        % plot bottom of cart
        set(h_carBase, 'XData', [rA(1) - p(6)/2, rA(1) + p(6)/2]);
        set(h_carBase, 'YData', [rA(2) - p(5)/2, rA(2) - p(5)/2]);
        
        % plot top of cart
        set(h_carTop, 'XData', [rA(1) - p(6)/2, rA(1) + p(6)/2]);
        set(h_carTop, 'YData', [rA(2) + p(5)/2, rA(2) + p(5)/2]);
        
        % plot link 1
        set(h_link1, 'XData', [rA(1) rB(1)]);
        set(h_link1, 'YData', [rA(2) rB(2)]);
        
        % plot link 2
        set(h_link2, 'XData', [rB(1) rC(1)]);
        set(h_link2, 'YData', [rB(2) rC(2)]);
        
        % plot link 3
        set(h_link3, 'XData', [rE(1) rD(1)]);
        set(h_link3, 'YData', [rE(2) rD(2)]);
        
        % plot ball
        set(ball, 'XData', [rEE(1, i) + ballX(1, i)*cosd(-theta(i)), rEE(1, i) + ballX(1, i)*cosd(-theta(i)) + r*sind(-theta(i))]);
        set(ball, 'YData', [rEE(2, i) - ballX(1, i)*sind(-theta(i)), rEE(2, i) - ballX(1, i)*sind(-theta(i)) + r*cosd(-theta(i))]);

        pause(0.1)
    end
end