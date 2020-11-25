function control_law = get_controller(zf, p_model, chosen_ctrl_str)
% get_controller: returns function handle that calculates control input
% according to the desired control law
%
% INPUTS
% zf: desired final state
% p_model: model parameters (may differ from 'real' system parameters)
% chosen_ctrl_str: the controller to be used
%
% OUTPUTS
% control_law: function handle to calculate control input at each step

switch chosen_ctrl_str
    
    % design lqr for feedback-linearized system
    case 'joint_space_fb_lin'
        nDof = 4;
        nState = nDof*2;
        
        A = [zeros(nDof) eye(nDof); zeros(nDof) zeros(nDof)];
        B = [zeros(nDof); eye(nDof)];
        C = eye(nState); % full state feedback for identity matrix, only measure (x,y)
        assert(rank(ctrb(A, B)) == length(A), 'System is not controllable\n'); % check controllability

        % Q matrix
        Q_fl = eye(nState)*100;
        Q_fl(2, 2) = 1000;

        % R matrix
        R_fl = [0.001 0.0   0.0   0.0;
                0.0   0.001 0.0   0.0;
                0.0   0.0   0.001 0.0;
                0.0   0.0   0.0   0.1];

        K_fl = lqr(A, B, Q_fl, R_fl); % lqr gains

        p_ctrl_fl.lqr_gain = K_fl;
        kr_q = -inv(C(1:4, 1:8)*inv(A - B*K_fl)*B);
        p_ctrl_fl.kr = kr_q; %-inv(C*inv(A - B*K_fl)*B); % scale the reference appropriately

    %     Kr = -inv(C*inv(A - B*K)*B);

        p_ctrl_fl.zf = zf; % desired final state

        control_law = @(t,z) control_law_feedback_linearization(t, z, p_model, p_ctrl_fl);

    % design lqr for feedback-linearized system in operational space
    case 'operational_space_fb_lin'
        nDof = 2;
        nState = nDof*2;
        
        A = [zeros(nDof) eye(nDof); zeros(nDof) zeros(nDof)];
        B = [zeros(nDof); eye(nDof)];
        C = eye(nState); % full state feedback for identity matrix, only measure (x,y)
        assert(rank(ctrb(A, B)) == length(A), 'System is not controllable\n'); % check controllability

        % Q matrix
        Q = eye(nState);
        Q(3, 3) = 10;
        Q(4, 4) = 10;
        
        % R matrix
%         R_fl = [0.001 0.0   0.0   0.0;
%                 0.0   0.001 0.0   0.0;
%                 0.0   0.0   0.001 0.0;
%                 0.0   0.0   0.0   0.1];
        R = [0.01 0.0;
             0.0  0.01];

        p_ctrl.lqr_gain = lqr(A, B, Q, R);
        p_ctrl.kr = eye(2); %-inv(C*inv(A-B*K)*B);
        p_ctrl.r = [0.5; 0.5]; % desired final state

        control_law = @(t,z) control_law_fb_lin_op_sp(t, z, p_model, p_ctrl);

    % design standard lqr towards desired final state
    case 'standard_lqr'
        % tangent linearization of dynamics about final desired position
        u_equi = Grav_arm(zf, p_model); % inputs at equilibrium = gravity - for non-zero set point
        [A_lin, B_lin] = linearize_dynamics(zf, u_equi, p_model);
        assert(rank(ctrb(A_lin, B_lin)) == length(A_lin), 'System is not controllable\n'); % check controllability

        C_lin = eye(length(A_lin));

        % Q matrix
        % TODO: Tune me!
        Q = zeros(8);
%         Q = C_lin'*C_lin*100;
        Q(1, 1) = 100;
        Q(2, 2) = 10; Q(3, 3) = 10; Q(4, 4) = 10;
        
        % R matrix
        % TODO: Tune me!
        R = eye(size(B_lin, 2))*1;

        K_lqr = lqr(A_lin, B_lin, Q, R);
        p_ctrl_lqr.lqr_gain = K_lqr;
        p_ctrl_lqr.zf = zf; % desired final state
        p_ctrl_lqr.kr = -(C_lin(1:4, :)/inv(A_lin - B_lin*K_lqr)*B_lin)*0; % set reference input to 0

        control_law = @(t, z) (control_law_standard_lqr(t, z, p_model, p_ctrl_lqr) + u_equi);

    otherwise
        error('Choose controller from available options')
end
    