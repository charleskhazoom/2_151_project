function simulate_arm()
% simulate_arm: calculates states and inputs over a time span and then
% plots the energy and kinematics of the system. Also animates the motion
% of the system

% INPUTS
% N/A

% OUTPUTS
% N/A
    
clear all; clc; %close all;

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
    
    %k = 1/2; % inertia of cylinder, kg-m^2
    k = 2/5; % inertia of ball, kg-m^2
   
    % Comment out ball mass, friction
%     m_ball = 0.1; % ball mass, kg
%     % seems like making mass larger will result in ball moving more
%     mu = 0.25; % coefficient of friction between ball and platform

%% Parameter vector (real system)
    p = [m_cart m1 m2 m3 h_cart l_cart l_1 l_2 l_3 g k]'; % parameters vector

%% Parameter vector (estimated system)

    p_estim = p; % estimate parameter vector
    mass_error = 1.05; % assume all masses are estimated incorrectly by some percentage
    p_estim(1:4) = p_estim(1:4)*mass_error;
    
%% Set up dynamic simulation    
    numStates = 10; % [x_cart, th_1, th_2, th_3, vx_cart, omega_1, omega_2, omega_3, ball_x, ball_v]'
    numInputs = 4; % [f_cart, t_joint1, t_joint2, t_joint3]'

    dt = 0.001; % timestep, sec

    tf = 6; % % final time, sec (change if not enough to complete task)
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step);
    
    p_cup_initial = [-0.8, 0.5]'; % initial position of ball (x, y), m
    q0 = eval(invKin_arm(p_cup_initial, p , [0, 0]')); % initial configuration
    z0 = [q0; 0; 0; 0; 0; 0; 0]; % initial state
    
    fprintf(['Initial State\nx_cart: ' num2str(q0(1)) ' m\n' ...
        'theta_1: ' num2str(q0(2)) ' rad\n' ...
        'theta_2: ' num2str(q0(3)) ' rad\n' ...
        'theta_3: ' num2str(q0(4)) ' rad\n']);

    p_cup_final = [0.5, 0.3]'; % need to specify orientation of last link in the world frame too!
    qf = eval(invKin_arm(p_cup_final, p, q0(2:3))); % final configuration
    zf = [qf; 0; 0; 0; 0; 0; 0]; % final state, consider if you can't overconstrain ball final state
    
    fprintf(['\nFinal State\nx_cart: ' num2str(qf(1)) ' m\n' ...
        'theta_1: ' num2str(qf(2)) ' rad\n' ...
        'theta_2: ' num2str(qf(3)) ' rad\n' ...
        'theta_3: ' num2str(qf(4)) ' rad\n']);
    
    z_out = zeros(numStates, num_step);
    z_out(:, 1) = z0;
    
    z_hat_out = zeros(numStates, num_step);
    z_hat_out(:,1) = z0*0.9; % initial observer states
    
    
    dz_out = zeros(numStates, num_step); % store rate of change of states at each time step
    dz_hat_out = zeros(numStates, num_step);% store rate of change of observer states at each time step

    u_out = zeros(numInputs, num_step); % store control input at each time step 
    
    z_int = zeros(2, num_step);% integral states for lqi controllers;
    z_int(:,1) = [0;0]; % initialize to zero; 
    
    ball_alongPlate = zeros(3, num_step); % ball position over time
    ball_alongPlate(:, 1) = [0; 0; 0];
    accel = zeros(2, num_step); % acceleration of platform
    %% Measurement matrix C
    
    Cob = zeros(5,numStates);
    % measure joint angles
    Cob(1,1) = 1;
    Cob(2,2) = 1;
    Cob(3,3) = 1;
    Cob(4,4) = 1;
    
    % measure ball position
    Cob(5,9) = 1;
    
    % I don't know how to implement a kalman filter here. I Cant' get the
    % kalman() to output a solution... Why?
    % noise covariance Wo
    Wo = zeros(5);
    Wo(1,1) = (1/12)*0.2;
    Wo(2,2) = (1/12);
    Wo(3,3) = (1/12);
    Wo(4,4) = (1/12);
    Wo(5,5) = 0.05^2;
    Wi=eye(4);
%     Wi(1,1) = 0.1^2;%;Wi(2,2) = 1^2;Wi(3,3) = 1^2;Wi(4,4) = 1^2;
    
    
    v = Wo*randn(5,num_step); % sensor noises
%% Choose control law
use_observer = 1;% 0: full state feedback. 1: observer feedback.
% also, if no observer is designed in the chosen control law, get
% controller returns obsv_dynamics as empty, and observer is ignored in
% integration loop

use_noise = 0;

%     ctrl_law_str = 'joint_space_fb_lin';
%     ctrl_law_str = 'joint_space_fb_lin_with_ball';
      ctrl_law_str = 'joint_space_fb_lin_with_ball_lqi';


%     ctrl_law_str = 'operational_space_fb_lin';
%     ctrl_law_str = 'standard_lqr';
    
    fprintf(['\nChosen control law: ' ctrl_law_str '\n\n']) 
    
    % outputs function handle to be used during Euler integration (for loop below)
    [control_law, obsv_dynamics,int_dynamics] = get_controller(zf, p_estim, ctrl_law_str,Cob,Wo,Wi);

% to design a new control law:
% 1) create function in external file which takes as argument
% t, z, u, p, p_ctrl (see example control_law_feedback_linearization.m)
% 2) edit get_controller, where the design parameters are created and
% passed to the anonymous function which is the output of control law.

%% Perform Euler Integration    
    for i = 1:num_step
        


        % choose which states to feed back
        if use_observer == 1 && ~isempty(obsv_dynamics)% feed back observer states
            z_ctrl = z_hat_out(:,i);
        else % full state feedback
            z_ctrl = z_out(:, i);
        end
        
        % compute integral states if any
        if ~isempty(int_dynamics)
            z_ctrl = [z_ctrl;z_int(:,i)]; % append integral states
            dz_int = int_dynamics(z_ctrl);
            
            %integral state (for lqi controllers)
            z_int(:,i+1) = dz_int*dt; % integrate error
        end
        
        
        
        % Compute Controls
        % get control input for this timestep
        u_out(:, i) = control_law(tspan(i), z_ctrl);
        % u_out(:, i) = zeros(4, 1);
        


        
        % calculate dz, change in state variables
        dz = dynamics(tspan(i), z_out(:, i), u_out(:, i), p);

        % Forward Euler to calculate next state
        z_out(:, i + 1) = z_out(:, i) + dz*dt;
        dz_out(:, i) = dz;
        
        % integrate observer dynamics
        if use_observer == 1 && ~isempty(obsv_dynamics)
            % measurements
            y = Cob*z_out(:,i) + v(:,i)*use_noise;

            dz_hat_out(:,i) = obsv_dynamics(y,z_hat_out(:,i), u_out(:, i));
            z_hat_out(:, i + 1) = z_hat_out(:, i) + dz_hat_out(:,i)*dt;
        end
        
        % Arent these lines redundant from 4 lines up?
        % Position update
%         z_out(1:4, i + 1) = z_out(1:4, i) + z_out(5:8, i + 1)*dt; % robot
%         z_out(9, i + 1) = z_out(9, i) + z_out(10, i + 1)*dt; % ball
        
        % Ball Simulation
        theta = z_out(2, i) + z_out(3, i) - 90/180*pi + z_out(4, i); % plate angle
        accel(:, i) = acceleration_endEffector([z_out(1:8, i); dz_out(5:8, i)], p);
        a_x_plate = accel(1, i)*cos(-theta) - accel(2, i)*sin(-theta); % rotate into plate frame
        
        a_b_plate = (1/(1 + k))*(g*sin(-theta) - a_x_plate); %
        
        % update ball kinematics
        ball_alongPlate(3, i + 1) = a_b_plate;
        ball_alongPlate(2, i + 1) = ball_alongPlate(2, i) + a_b_plate*dt;
        ball_alongPlate(1, i + 1) = ball_alongPlate(1, i) + ball_alongPlate(2, i + 1)*dt;

    end
    
    % final states
    z_out = z_out(:, 1:end - 1);
    ball_alongPlate = ball_alongPlate(:, 1:end - 1);
    
%% Look at Results
    plot_obsv = use_observer && ~isempty(obsv_dynamics);
    make_plots(tspan, z_out, u_out, dz_out,z_hat_out,dz_hat_out, ball_alongPlate, accel, p,plot_obsv);

    rE = zeros(2, length(tspan)); % end effector position
    vE = zeros(2, length(tspan)); % end effector velocity
    
    for i = 1:length(tspan)
        rE(:, i) = position_endEffector(z_out(:, i), p);
        vE(:, i) = velocity_endEffector(z_out(:, i), p);
    end

%% Animate Solution
    figure(7); clf;
    theta = z_out(2, :) + z_out(3, :) - 90/180*pi + z_out(4, :); % plate angle
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
    
    ballPlot = plot([0], [0], 'LineWidth', 2); % ball position
    r = 0.1; % ball radius, m (arbitrary)    
    
    xlabel('x'); ylabel('y');
    h_title = title('t = 0.0 s');
    
    axis equal
    axis([-1 1 -1 1]);

    % Step through and update animation
    for i = 1:length(tspan)
      
        % skip eac 50th frame.
        if mod(i,50)
            continue;
        end
        
        t = tspan(i); % time
        z = x(:, i); % state
        keypoints = keypoints_arm(z, p); % defining parts of arm

        rA = real(keypoints(:, 1)); % center of mass (which is at 0) of cart
        rB = real(keypoints(:, 2)); % where link 1 and 2 meet
        rC = real(keypoints(:, 3)); % where link 2 and 3 meet
        rD = real(keypoints(:, 4)); % right side of plate(3)
        rE = real(keypoints(:, 5)); % left side of plate(3)

        set(h_title,'String',  sprintf('t = %.2f', t) ); % update title (time)
        
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
        set(h_link1, 'XData', [rA(1), rB(1)]);
        set(h_link1, 'YData', [rA(2), rB(2)]);
        
        % plot link 2
        set(h_link2, 'XData', [rB(1), rC(1)]);
        set(h_link2, 'YData', [rB(2), rC(2)]);
        
        % plot link 3
        set(h_link3, 'XData', [rE(1), rD(1)]);
        set(h_link3, 'YData', [rE(2), rD(2)]);
        
        % plot ball
        % ball = ballX(1, i); % check if it works with manual simulation of ball
         ball = z(9); % ball position
        xcenter  = mean([rEE(1, i) + ball*cosd(-theta(i)), rEE(1, i) + ball*cosd(-theta(i)) + r*sind(-theta(i))]);
        ycenter  = mean([rEE(2, i) - ball*sind(-theta(i)), rEE(2, i) - ball*sind(-theta(i)) + r*cosd(-theta(i))]);
        
%         set(ballPlot,'XData',[rEE(1, i) + ball*cosd(-theta(i)), rEE(1, i) + ball*cosd(-theta(i)) + r*sind(-theta(i))])
%         set(ballPlot,'YData',[rEE(2, i) - ball*sind(-theta(i)), rEE(2, i) - ball*sind(-theta(i)) + r*cosd(-theta(i))])

        set(ballPlot, 'XData', real(xcenter),'Marker','o','MarkerFaceColor','r','color','r','MarkerSize',10);
        set(ballPlot, 'YData', real(ycenter));%,'Marker','o','MarkerFaceColor','r','color','r');

        pause(0.1) % wait, draw next frame
    end
end