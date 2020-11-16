function control_law = get_controller(zf,p_model,chosen_crtl_str)

% p_model are model parameters, which can be different from the real system
% parameters

switch chosen_crtl_str
    case 'joint_space_fb_lin'
        %% design lqr for feedback-linearized system
    A_fl = [zeros(4,4) eye(4,4);zeros(4,4)  zeros(4,4)];
    B_fl = [zeros(4,4);eye(4,4)];
    C_fl = eye(8,8); % full state feedback for identity matrix
    rank(ctrb(A_fl,B_fl)) % full rank
    
    Q_fl = [eye(8,8)]*100;
    Q_fl(2,2) = 1000;
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
    
    control_law = @(t,z) control_law_feedback_linearization(t, z, p_model,p_ctrl_fl);
    
    case 'operational_space_fb_lin'
        nDof=2;
        nState = nDof*2;
        
    A = [zeros(nDof,nDof) eye(nDof,nDof);zeros(nDof,nDof)  zeros(nDof,nDof)];
    B= [zeros(nDof,nDof);eye(nDof,nDof)];
    C= eye(nState,nState); % full state feedback for identity matrix, % only measure x,y
    
    rank(ctrb(A,B)) % full rank
    
    Q= [eye(nState,nState)];
    Q(3,3)=10;
    Q(4,4)=10;
    %     Q_fl(2,2) = 1000;
%     R= [0.001 0 0 0 ;...
%             0 0.001 0 0;...
%             0 0 0.001 0;...
%             0 0 0     0.1];
    R = [0.01 0;...
         0 0.01];


    p_ctrl.lqr_gain = lqr(A,B,Q,R);
            
                
    p_ctrl.kr = eye(2);%-inv(C*inv(A-B*K)*B);
    
    p_ctrl.r = [0.5;0.5]; % desired final state

    
    control_law = @(t,z) control_law_fb_lin_op_sp(t, z, p_model,p_ctrl);
    
    case 'standard_lqr'
    % tangent linearization of dynamics about final desired position
    u_equi = Grav_arm(zf,p_model); % inputs at equilibrium = gravity
    [A_lin, B_lin] =  linearize_dynamics(zf,u_equi,p_model);
    
    C_lin = eye(size(A_lin, 1));
    
    Q = C_lin'*C_lin*1000; % TODO: Tune me!
    R = eye(size(B_lin, 2))*1; % TODO: Tune me!
    
    K_lqr = lqr(A_lin, B_lin, Q, R);
    p_ctrl_lqr.lqr_gain = K_lqr;
    p_ctrl_lqr.zf = zf; % desired final state
    p_ctrl_lqr.kr = -inv(C_lin(1:4, :)*inv(A_lin - B_lin*K_lqr)*B_lin)*0;
    
    control_law = @(t, z) (control_law_standard_lqr(t, z, p_model, p_ctrl_lqr)+u_equi);
    
    otherwise
        error('choose controller among available options')
end
    