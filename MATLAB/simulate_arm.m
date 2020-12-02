function simulate_arm()
% simulate_arm: calculates states and inputs over a time span and then
% plots the energy and kinematics of the system. Also animates the motion
% of the system

% INPUTS
% N/A

% OUTPUTS
% N/A
    
clear; clc; %close all;

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
    
    pos_ball0 = [-0.8, 0.5]'; % initial position of ball (x, y), m
    q0 = eval(invKin_arm(pos_ball0, p , [0, 0]')); % initial configuration
    z0 = [q0; 0; 0; 0; 0; 0.02; 0]; % initial state
    
    fprintf(['Initial State\nx_cart: ' num2str(q0(1)) ' m\n' ...
        'theta_1: ' num2str(q0(2)) ' rad\n' ...
        'theta_2: ' num2str(q0(3)) ' rad\n' ...
        'theta_3: ' num2str(q0(4)) ' rad\n']);

    pos_ballf = [0.5, 0.3]'; % need to specify orientation of last link in the world frame too!
    qf = eval(invKin_arm(pos_ballf, p, q0(2:3))); % final configuration
    zf = [qf; 0; 0; 0; 0; 0; 0]; % final state, consider if you can't overconstrain ball final state
    
    fprintf(['\nFinal State\nx_cart: ' num2str(qf(1)) ' m\n' ...
        'theta_1: ' num2str(qf(2)) ' rad\n' ...
        'theta_2: ' num2str(qf(3)) ' rad\n' ...
        'theta_3: ' num2str(qf(4)) ' rad\n']);
    
    z_out = zeros(numStates, num_step);
    z_out(:, 1) = z0;
    
    z_hat_out = zeros(numStates, num_step);
    z_hat_out(:, 1) = z0*0.95; % initial observer states
    
    dz_out = zeros(numStates, num_step); % store rate of change of states at each time step
    dz_hat_out = zeros(numStates, num_step); % store rate of change of observer states at each time step

    u_out = zeros(numInputs, num_step); % store control input at each time step 
    
    z_int = zeros(2, num_step); % integral states for lqi controllers;
    z_int(:, 1) = [0; 0]; % initialize to zero; 
    
    ball_alongPlate = zeros(3, num_step); % ball position over time
    ball_alongPlate(:, 1) = [0; 0; 0];
    accel = zeros(2, num_step); % acceleration of platform
    
%% Measurement matrix C
    % measure cart position, joint angles, ball position
    Cob = zeros(5, numStates);
    Cob(1, 1) = 1;
    Cob(2, 2) = 1;
    Cob(3, 3) = 1;
    Cob(4, 4) = 1;
    Cob(5, 9) = 1;
    
    % I don't know how to implement a kalman filter here. I Cant' get the
    % kalman() to output a solution... Why?
    % noise covariance Wo
    Wo = zeros(size(Cob, 1));
    Wo(1, 1) = (1/12)*0.2;
    Wo(2, 2) = (1/12);
    Wo(3, 3) = (1/12);
    Wo(4, 4) = (1/12);
    Wo(end, end) = (0.05)^2;
    
    Wi = eye(4)*0.1^2;
%     Wi(1, 1) = 0.1^2; Wi(2, 2) = 1^2; Wi(3, 3) = 1^2; Wi(4, 4) = 1^2;
    
    v = Wo*randn(size(Cob, 1), num_step); % sensor noises
    
%% Choose control law
    use_observer = 1; % 0: full state feedback. 1: observer feedback.
    % also, if no observer is designed in the chosen control law, get
    % controller returns obsv_dynamics as empty, and observer is ignored in
    % integration loop

    use_noise = 0; % 0: no noise. 1: Gaussian noise.
    
    use_mass_error = 1; % 0: no mass error. 1: mass error.

%     ctrl_law_str = 'joint_space_fb_lin';
%     ctrl_law_str = 'joint_space_fb_lin_with_ball';
%     ctrl_law_str = 'joint_space_fb_lin_with_ball_lqi';

%     ctrl_law_str = 'operational_space_fb_lin';
    ctrl_law_str = 'standard_lqr';
    
    fprintf(['\nChosen control law: ' ctrl_law_str '\n']) 
    fprintf(['Using observer: ' num2str(use_observer) '\n'])
    fprintf(['Using noise: ' num2str(use_noise) '\n'])
    
    if use_mass_error
        params = p_estim;
    else
        params = p;
    end        
    
    % outputs function handle to be used during Euler integration (for loop below)
    [control_law, obsv_dynamics, int_dynamics] = get_controller(zf, params, ctrl_law_str, Cob, Wo, Wi);

% to design a new control law:
% 1) create function in external file which takes as argument
% t, z, u, p, p_ctrl (see example control_law_feedback_linearization.m)
% 2) edit get_controller, where the design parameters are created and
% passed to the anonymous function which is the output of control law.

%% Perform Euler Integration    
    for i = 1:num_step

        % choose which states to feed back
        if use_observer && ~isempty(obsv_dynamics) % feed back observer states
            z_ctrl = z_hat_out(:, i);
        else % full state feedback
            z_ctrl = z_out(:, i);
        end
        
        % compute integral states (if any)
        if ~isempty(int_dynamics)
            z_ctrl = [z_ctrl; z_int(:, i)]; % append integral states
            dz_int = int_dynamics(z_ctrl);
            
            %integral state (for lqi controllers)
            z_int(:, i + 1) = dz_int*dt; % integrate error
        end
        
        % Compute Controls
        % get control input for this timestep
        u_out(:, i) = control_law(tspan(i), z_ctrl);
        
        % calculate dz, change in state variables
        dz = dynamics(tspan(i), z_out(:, i), u_out(:, i), p);

        % Forward Euler to calculate next state
        z_out(:, i + 1) = z_out(:, i) + dz*dt;
        dz_out(:, i) = dz;
        
        % integrate observer dynamics
        if use_observer && ~isempty(obsv_dynamics)
            % measurements
            y = Cob*z_out(:, i) + v(:, i)*use_noise;

            dz_hat_out(:, i) = obsv_dynamics(y, z_hat_out(:, i), u_out(:, i));
            z_hat_out(:, i + 1) = z_hat_out(:, i) + dz_hat_out(:, i)*dt;
        end
        
        % Ball Simulation
        theta = z_out(2, i) + z_out(3, i) - pi/2 + z_out(4, i); % plate angle
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
    make_plots(tspan, z_out, u_out, dz_out, z_hat_out, dz_hat_out, ball_alongPlate, accel, p, plot_obsv, zf);

    rE = zeros(2, length(tspan)); % end effector position
    vE = zeros(2, length(tspan)); % end effector velocity
    
    for i = 1:length(tspan)
        rE(:, i) = position_endEffector(z_out(:, i), p);
        vE(:, i) = velocity_endEffector(z_out(:, i), p);
    end

%% Animate Solution
    % option to save video
    save_video = 0;
    % get video title
    videoName = [ctrl_law_str '_observer' num2str(use_observer) '_mass_error' num2str(use_mass_error)];
    
    figure(7); clf;
    theta = z_out(2, :) + z_out(3, :) - pi/2 + z_out(4, :); % plate angle
    hold on
    keep_frames = 0;
    animateSol(tspan, z_out, p, ball_alongPlate, rE, theta, pos_ball0, pos_ballf, keep_frames, videoName, save_video);
    
end