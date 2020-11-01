function simulate_coffeeArm()
    clear; clc
    %% Define fixed paramters
    m1 =50 ;         
    m2 =5; 
    m3 =5;            
    m4 =5;
    h_1 = 6*0.0254;      
    l_1 = 12*0.0254;
    l_2 = 12*0.0254;
    l_3 = 12*0.0254;
    l_4 = 12*0.0254;
    g = 9.81;    

    %% Parameter vector (real system)
    p   = [m1 m2 m3 m4 h_1 l_1 l_2 l_3 l_4 g]';        % parameters
    %% Parameter vector (estimated system)

    p_estim = p; % e
    p_estim(1:4)=p_estim(1:4)*1; % assume all masses are overestimated by some percentage
    %% Simulation Parameters Set 2 -- Operational Space Control
    
    p_traj.omega = 3;
    p_traj.x_0   = 0;
    p_traj.y_0   = -.125;
    p_traj.r     = 0.025;
    
    
    %% Perform Dynamic simulation
    dt = 0.001;
    tf = 2; %May have to change if 10 second not enough to complete task
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 
    
    p_cup_initial = [0.1,0.85]';
    q0 = invKin_arm(p_cup_initial,p,[0,0,0,0]');
    z0 = [q0;0;0;0;0];
    
    p_cup_final = [3,0.2]'; % need to specify orientation of last link in the world frame too!
    qf = eval(invKin_arm(p_cup_final,p,q0));
    zf = [qf;0;0;0;0];
    
    z_out = zeros(8,num_step);
    z_out(:,1) = z0;
    
    %% design lqr for feedback-linearized system
    A_fl = [zeros(4,4) eye(4,4);zeros(4,4)  zeros(4,4)];
    B_fl = [zeros(4,4);eye(4,4)];
    C_fl = eye(8,8); % full state feedback for identity matrix
    rank(ctrb(A_fl,B_fl)) % full rank
    
    Q_fl = [eye(8,8)]*1000;
    R_fl = [0.001 0 0 0 ;...
            0 0.001 0 0;...
            0 0 0.001 0;...
            0 0 0     0.1];
    K_fl = lqr(A_fl,B_fl,Q_fl,R_fl);
%     MassMatrix = A_
    p_ctrl_fl.lqr_gain = K_fl;
    kr_q = -inv(C_fl(1:4,1:8)*inv(A_fl-B_fl*K_fl)*B_fl);
    
    p_ctrl_fl.kr = [kr_q];%-inv(C_fl*inv(A_fl-B_fl*K_fl)*B_fl); % scale the reference appropriately
%     Kr = -inv(C*inv(A-B*K)*B);
    p_ctrl_fl.zf = zf; % desired final state
    %%
    
    for i=1:num_step-1
        u_out(:,i) = control_law_feedback_linearization(tspan(i),z_out(:,i),p_estim,p_ctrl_fl);
        dz = dynamics(tspan(i), z_out(:,i),u_out(:,i), p, p_traj);
        %z_out(3,i) = joint_limit_constraint(z_out(:,i),p);
        z_out(:,i+1) = z_out(:,i) + dz*dt;
        
        % Position update
        z_out(1:4,i+1) = z_out(1:4,i) + z_out(5:8,i+1)*dt;
    end
    
    %% Compute Energy (shouldnt need changing)
    E = energy_coffeeArm(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)');
    
    %% Compute foot position over time (havent look in depth)
    rE = zeros(2,length(tspan));
    vE = zeros(2,length(tspan));
    for i = 1:length(tspan)
        rE(:,i) = position_endEffector(z_out(:,i),p);
        vE(:,i) = velocity_endEffector(z_out(:,i),p);
    end
    
    figure(2); clf;
    plot(tspan,rE(1,:),'r','LineWidth',2)
    hold on
    plot(tspan,p_traj.x_0 + p_traj.r * cos(p_traj.omega*tspan) ,'r--');
    plot(tspan,rE(2,:),'b','LineWidth',2)
    plot(tspan,p_traj.y_0 + p_traj.r * sin(p_traj.omega*tspan) ,'b--');
    
    
    xlabel('Time (s)'); ylabel('Position (m)'); legend({'x','x_d','y','y_d'});

    figure(3); clf;
    plot(tspan,vE(1,:),'r','LineWidth',2)
    hold on
    plot(tspan,vE(2,:),'b','LineWidth',2)
    
    xlabel('Time (s)'); ylabel('Velocity (m)'); legend({'vel_x','vel_y'});
    
    figure(4)
    plot(tspan,z_out(1:2,:)*180/pi)
    legend('q1','q2');
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    
    figure(5)
    plot(tspan,z_out(3:4,:)*180/pi)
    legend('q1dot','q2dot');
    xlabel('Time (s)');
    ylabel('Angular Velocity (deg/sec)');
    
    %% Animate Solution
    figure(6); clf;
    hold on
   
  
    % Target traj
    TH = 0:.1:2*pi;
    plot( p_traj.x_0 + p_traj.r * cos(TH), ...
          p_traj.y_0 + p_traj.r * sin(TH),'k--'); 
    
    animateSol(tspan, z_out,p);
end

function u = control_law(t, z, p, p_traj)
u = [0; 0; 0; 0];
%     % Controller gains, Update as necessary for Problem 1
%     K_x = 150.; % Spring stiffness X
%     K_y = 150.; % Spring stiffness Y
%     D_x = 10.;  % Damping X
%     D_y = 10.;  % Damping Y
% 
%     %will be useful if we make a reference trajectory
% %     % Desired position of foot is a circle
% %     omega_swing = p_traj.omega;
% %     rEd = [p_traj.x_0 p_traj.y_0 0]' + ...
% %             p_traj.r*[cos(omega_swing*t) sin(omega_swing*t) 0]';       
%     
%     % Actual position and velocity 
%     rE = position_endEffector(z,p);
%     vE = velocity_endEffector(z,p);
%     
%     % Quesiton 1.1
%     J  = jacobian_foot(z,p);
%     Mass = A_coffeeArm(z,p);
%     V = Corr_arm(z,p); 
%     jdot = jacobian_dot_endEffector(z,p);
%     G = Grav_arm(z,p);
%     
%     A = inv(J*inv(Mass)*J');
%     mu = A*J*inv(Mass)*V-A*jdot*z(3:4);
%     % 1.4 Improper estimation of mu
%     %random= 1+((rand(1)*2-2)*.2); %+- 20%
%     %mu = random*(A*J*inv(Mass)*V-A*jdot*z(3:4));
%     
%     rho = A*J*inv(Mass)*G;
%     
%     f  = [aEd(1)+ K_x * (rEd(1) - rE(1) ) + D_x * (vEd(1) - vE(1) ) ;
%           aEd(2)+ K_y * (rEd(2) - rE(2) ) + D_y * (vEd(2) - vE(2) ) ];
%    
%     tau = J' *(A*f + mu + rho);
%     
%     % Compute virtual foce for Question 1.4 and 1.5
%     %f  = [K_x * (rEd(1) - rE(1) ) + D_x * ( - vE(1) ) ;
%     %      K_y * (rEd(2) - rE(2) ) + D_y * ( - vE(2) ) ];
%     
%     %% Task-space compensation and feed forward for Question 1.8
% 
%     % Map to joint torques  
%     %tau = J' * f;
end
function u = control_law_feedback_linearization(t, z, p,p_ctrl_fl)
      u = [0; 0; 0; 0];
%   unpack controller parameters
    
    K_lqr = p_ctrl_fl.lqr_gain;
    kr = p_ctrl_fl.kr; % scale the reference appropriately
    zf = p_ctrl_fl.zf; % desired final state
%     % Actual position and velocity 
%     
%     % Quesiton 1.1
%     J  = jacobian_foot(z,p);
    M = A_coffeeArm(z,p);
    V = Corr_arm(z,p); 
%     jdot = jacobian_dot_endEffector(z,p);
    G = Grav_arm(z,p);
    
    r  = zf(1:4);% desired new equilibrium point
%     kr = % gain to 
    w = kr*r-K_lqr*z; % do lqr here 
    u = V+G + M*w;

end

function dz = dynamics(t,z,u,p,p_traj)
    % Get mass matrix
    A = A_coffeeArm(z,p);
    
    % Compute Controls % HERE WE CAN CHANGE THE CONTROL LAW TO BE ANYTHING
%     u = control_law(t,z,p,p_traj);
%     u = control_law_feedback_linearization(t,z,p_estim,p_ctrl);
    
    % Get b = Q - V(q,qd) - G(q)
    b = b_coffeeArm(z,u,p);
    
    % Solve for qdd.
    qdd = A\(b);
    dz = 0*z;
    
    % Form dz
    dz(1:4) = z(5:8);
    dz(5:8) = qdd;
end

function qdot = discrete_impact_contact(z,p, rest_coeff, fric_coeff, yC)
    rE = position_endEffector(z,p);
    C = rE(2) - yC;
    vE = velocity_endEffector(z,p);
    Cdot = vE(2);
    J = jacobian_foot(z,p);
    Mass = A_coffeeArm(z,p);
    A = inv(J*inv(Mass)*J');
    
    if C<0 & Cdot<0
        Fcz = [0 0; A(2,:)]*(-rest_coeff*Cdot - [0 0; J(2,:)]*z(3:4));
        qdot = z(3:4) + inv(Mass)*J(2,:)'.*Fcz; %matrix size?
        Fcx = [A(1,:); 0 0]*(0-[J(1,:); 0 0]*z(3:4)); 
        if abs(Fcx(1)) > abs(fric_coeff*Fcz(2))
            Fcx(1) = -fric_coeff*Fcz(2);
        end
        qdot = qdot + inv(Mass)*J(1,:)'.*Fcx;
    else
        qdot = z(3:4);
    end
end

function qdot = joint_limit_constraint(z,p)
    C = z(1);
    Cdot = z(3);
    %Mass = A_coffeeArm(z,p);
    
    if C<-50/360*2*pi && Cdot<0
        qdot = z(3) + (-.1*Cdot-.1*C);
    else
        qdot = z(3);
    end
end

function animateSol(tspan, x, p)
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
        keypoints = keypoints_coffeeArm(z,p);

        rA = keypoints(:,1); % center of mass (which is at 0) of cart
        rB = keypoints(:,2); % where link 1 and 2 meet
        rC = keypoints(:,3); % where link 2 and 3 meet
        rD = keypoints(:,4); % where link 3 end

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
        
        set(h_link3,'XData',[rC(1) rD(1)]);
        set(h_link3,'YData',[rC(2) rD(2)]);

        pause(.01)
    end
end


function [Jx, Ju] =  linearize_dynamics(z_equi,u_equi,params)

         % Get mass matrix
    A = A_coffeeArm(z_equi,p);
    

%     u = control_law_feedback_linearization(t,z,p_estim,p_ctrl);
    
    % Get b = Q - V(q,qd) - G(q)
    b = b_coffeeArm(z_equi,u_equi,p);
    
    % Solve for qdd.
    qdd = A\(b);
    dz = 0*z;
    
    % Form dz
    dz(1:4) = z_equi(5:8);
    dz(5:8) = qdd;
%     side =1;

%     tspan = [0 1];
%     [~, zout] = slip_simulation_ode45_v2(z_equi,u_equi,params,side,tspan);
    
    
    zout0 = zout(:,end);

    %statef0
    ep = 1e-6;
    delStateMat = [0 0 1 0 0 0;0 0 0 1 0 0; 0 0 0 0 1 0]';
    Jx = zeros(3,3);
    Ju = zeros(3,3);

    for i=1:3
        z_dev = z_equi+ep*delStateMat(:,i);
        
        [~, zout] = slip_simulation_ode45_v2(z_dev,u_equi,params,side,tspan);  
        zoutf = zout(:,end);
        
        Jx(:,i) = (zoutf(3:5) - zout0(3:5))/ep;
    end


    del_u_mat= eye(3);
    del_u_mat = [del_u_mat;0 0 -1];
%     del_u_mat(end,end)= -1; % change in spring stiffness during decompression phase = -1*change in spring stiffness at compression phase
    
%     ep_u = [ep;ep;ep*1e4;ep*1e4];
    ep_u = [ep;ep;ep*1e4];
    
    for i=1:3
        
        u_dev = u_equi + del_u_mat(:,i)*ep_u(i);

%         dynParams.tdParams = [newParams(1) newParams(2)];
%         dynParams.heightThreshold = robotLegLength*cos(dynParams.tdParams(1));
%         dynParams.L0 = sqrt(robotLegLength^2+robotHipDisp^2 +...
%             2*robotLegLength*robotHipDisp*sin(dynParams.tdParams(1))*sin(dynParams.tdParams(2)));
%         dynParams.k1 = (dynParams.k0 + newParams(3))*1e4;
%         dynParams.k2 = (dynParams.k0 - newParams(3))*1e4;
        [~, zout] = slip_simulation_ode45_v2(z_equi,u_dev,params,side,tspan);
        zoutf = zout(:,end);
        
%         [T_out STATE_out FOOT_out tf statef EN_out stanceTime tdPos] = simulatePeriod(0,stateInit,dynParams);

        Ju(:,i) = (zoutf(3:5) - zout0(3:5))/ep_u(i);
    end

end