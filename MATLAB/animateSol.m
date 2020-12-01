function animateSol(tspan, x, p, ballX, rEE, theta, start_pos, final_pos, keep_frames)
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
% keep_frames: boolean about whether holding plots

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
    axis([-1 1 -0.1 1]);
     
    seg1x = [];
    seg1y = [];
    seg2x = [];
    seg2y = [];

    seg3x = [];
    seg3y = [];

    seg4x = [];
    seg4y = [];
    seg5x = [];
    seg5y = [];

    seg6x = [];
    seg6y = [];
    seg7x = [];
    seg7y = [];
    ball_traj = [];
    col = [0.3 0.3 0.3];
     
    % Step through and update animation
    k = 1;
    pause(1) % wait 1 second
    for i = 1:length(tspan)
        
        % skip eac 50th frame.
        if mod(i, 50)
            continue;
        end
        
        k = k + 1;
        t = tspan(i); % time
        z = x(:, i); % state
        
        if keep_frames    
            for j = 1:length(seg1x)/2
    %             plot(seg1x(2*j - 1:2*j), seg1y(2*j - 1:2*j), 'color', col); hold on
    %             plot(seg2x(2*j - 1:2*j), seg2y(2*j - 1:2*j), 'color', col); hold on
    %             plot(seg3x(2*j - 1:2*j), seg3y(2*j - 1:2*j), 'color', col); hold on
    %             plot(seg4x(2*j - 1:2*j), seg4y(2*j - 1:2*j), 'color', col); hold on
                plot(seg5x(2*j - 1:2*j), seg5y(2*j - 1:2*j), 'color', col); hold on
                plot(seg6x(2*j - 1:2*j), seg6y(2*j - 1:2*j), 'color', col); hold on
                plot(seg7x(2*j - 1:2*j), seg7y(2*j - 1:2*j), 'color', col); hold on
            end
        end 

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
        ball_traj = [ball_traj [xcenter; ycenter]];
%         set(ballPlot,'XData',[rEE(1, i) + ball*cosd(-theta(i)), rEE(1, i) + ball*cosd(-theta(i)) + r*sind(-theta(i))])
%         set(ballPlot,'YData',[rEE(2, i) - ball*sind(-theta(i)), rEE(2, i) - ball*sind(-theta(i)) + r*cosd(-theta(i))])

        set(ballPlot, 'XData', real(xcenter), 'Marker', 'o', 'MarkerFaceColor', 'r', 'color', 'r', 'MarkerSize', 12);
        set(ballPlot, 'YData', real(ycenter));%, 'Marker', 'o', 'MarkerFaceColor', 'r', 'color', 'r');
        
        if keep_frames
            if ~mod(k, 1) && tspan(i) < 0.5 || ~mod(k, 8) && tspan(i) >= 0.5 % keep some frames(i == 1)||(i == 20)||(i == 30)
                seg1x = [seg1x h_carLSide.XData];
                seg1y = [seg1y h_carLSide.YData];
                
                seg2x = [seg2x h_carRSide.XData];
                seg2y = [seg2y h_carRSide.YData];
                
                seg3x = [seg3x h_carBase.XData];
                seg3y = [seg3y h_carBase.YData];
                
                seg4x = [seg4x h_carTop.XData];
                seg4y = [seg4y h_carTop.YData];
                
                seg5x = [seg5x h_link1.XData];
                seg5y = [seg5y h_link1.YData];
                
                seg6x = [seg6x h_link2.XData];
                seg6y = [seg6y h_link2.YData];
                seg7x = [seg7x h_link3.XData];
                seg7y = [seg7y h_link3.YData];
            end
            
            plot(ball_traj(1, :), ball_traj(2, :), '--k', 'linewidth', 1.4);
        end
        
        pause(0.1) % wait, draw next frame
    end
end