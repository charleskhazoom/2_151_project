function [control_law, observer_dynamics, int_dynamics] = get_controller(zf, p_model, chosen_ctrl_str, C_ob, Wo, Wi)
% get_controller: returns function handles that calculates control input
% according to the desired control law, state estimate dynamics, and
% integrator dynamics
%
% INPUTS
% zf: desired final state
% p_model: model parameters (may differ from 'real' system parameters)
% chosen_ctrl_str: the controller to be used
% Cob: C matrix for an oersver - we can measure the 3 joint angles and the
% ball position on the plate
% Wo: output noise (standard Kalman formulation)
% Wi: input noise (standard Kalman formulation)
%
% OUTPUTS
% control_law: function handle to calculate control input at each step
% observer_dynamics: function handle to calculate observer dynamics, to be integrated in the main loop.
% int_dyanamics: integral dyanmics for integrators (for lqi for instance)


% initialize as empty. Will be overwritten by a fcn 
% handle only if an observer is designed/integral states are added (lqi).
observer_dynamics = [];
int_dynamics = [];

switch chosen_ctrl_str
    
    case 'joint_space_fb_lin'
        % design lqr for feedback-linearized system
        nDof = 4;
        nStates = nDof*2;
        
        zf = zf(1:nStates); % ignore ball states
        
        A = [zeros(nDof) eye(nDof); zeros(nDof) zeros(nDof)];
        B = [zeros(nDof); eye(nDof)];
        C = eye(nStates); % full state feedback for identity matrix, only measure (x,y)
        assert(rank(ctrb(A, B)) == length(A), 'System is not controllable\n'); % check controllability

        % Q matrix
        Q_fl = eye(nStates)*100;
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
        
    case 'joint_space_fb_lin_with_ball'
    % design lqr for feedback-linearized system
        nDof = 5;
        nInputs = 4;
        
        nStates = nDof*2;
        
        % feedback linearized dyanmics of arm only
        Afl = [zeros(nDof - 1) eye(nDof - 1); zeros(nDof - 1) zeros(nDof - 1)];
        Bfl = [zeros(nInputs); eye(nInputs)];
        
        % linearize ball dynamics wrt input qdd and state
        u_equi = [0; 0; 0; 0]; % inputs at equilibrium = no joint acceleration
        [Aball, Bball] = linearize_ball_dynamics(zf, u_equi, p_model);
        
        A = [Afl zeros(8, 2); Aball];
        B = [Bfl; Bball];
        
%         C = eye(nState); % full state feedback for identity matrix, only measure (x,y)
        assert(rank(ctrb(A, B)) == length(A), 'System is not controllable\n'); % check controllability

        % Q matrix
        Q = zeros(nStates);%*(1/deg2rad(10))^2;
        
        Q(1, 1) = 1/(0.1)^2;
        Q(2, 2) = 1/(deg2rad(10))^2;
        Q(3, 3) = 1/(deg2rad(10))^2;
        Q(4, 4) = 1/(deg2rad(10))^2;
        
        Q(5, 5) = 1/(2)^2; % forward cart velocity 
        Q(6, 6) = 1/(50)^2; % 30 rad/s ~= 300 rpm
        Q(7, 7) = 1/(50)^2;
        Q(8, 8) = 1/(50)^2;

        
        Q(9, 9) = 1/(0.1)^2;
        Q(10, 10) = 1/(1)^2;

%         Q(2, 2) = 1000;

        % R matrix
        R = zeros(nInputs);

        R(1, 1) = 1/(0.3)^2;
        R(2, 2) = 1/(0.5)^2;
        R(3, 3) = 1/(3)^2;
        R(4, 4) = 1/(10)^2;


        K = lqr(A, B, Q, R); % lqr gains

        p_ctrl_fl.lqr_gain = K;
        kr_q = zeros(4);%-inv(C([1:4], 1:10)*inv(A - B*K)*B);
        p_ctrl_fl.kr = kr_q; %-inv(C*inv(A - B*K_fl)*B); % scale the reference appropriately

    %     Kr = -inv(C*inv(A - B*K)*B);

        p_ctrl_fl.zf = zf; % desired final state

        control_law = @(t,z) control_law_feedback_linearization_with_ball(t, z, p_model, p_ctrl_fl);
        
    %--------------------------- observer design ---------------------------- %
        obs_poles = eig(A - B*K)*3; % arbitrary decision - poles 3x faster than system
     
        % pole placement (replace by Kalman later)
        L = place(A', C_ob', obs_poles)';
        p_obsv.L = L;
        p_obsv.A = A;
        p_obsv.B = B;
        p_obsv.C = C_ob;
        p_obsv.zf = zf;
        observer_dynamics = @(y, x_hat, u) obsv_feedback_linearization_with_ball(y, x_hat, u, p_model, p_obsv);

    case 'joint_space_fb_lin_with_ball_lqi'
        % design lqr for feedback-linearized system
        nDof = 5;
        nInputs = 4;
        
        nStates = nDof*2;
        
        % feedback linearized dyanmics of arm only
        Afl = [zeros(nDof - 1) eye(nDof - 1); zeros(nDof - 1) zeros(nDof - 1)];
        Bfl = [zeros(nInputs); eye(nInputs)];
        
        % linearize ball dynamics wrt input qdd and state
        u_equi = [0; 0; 0; 0]; % inputs at equilibrium = no joint acceleration
        [Aball, Bball] = linearize_ball_dynamics(zf, u_equi, p_model);
        
        A = [Afl zeros(8, 2); Aball];
        B = [Bfl; Bball];
        
        assert(rank(ctrb(A, B)) == length(A), 'System is not controllable\n'); % check controllability

        % Q matrix
        Q = eye(nStates + 2); % 2 additional integral states
        
        Q(1, 1) = 1/(0.1^2);
        Q(2, 2) = 1/(deg2rad(90))^2;
        Q(3, 3) = 1/(deg2rad(90))^2;
        Q(4, 4) = 1/(deg2rad(90))^2;
        
        Q(5, 5) = 1/(10)^2; % forward cart velocity 
        Q(6, 6) = 1/(80)^2; % 30 rad/s ~= 300 rpm
        Q(7, 7) = 1/(80)^2;
        Q(8, 8) = 1/(80)^2;
        
        Q(9, 9) = 1/(0.2)^2;
        Q(10, 10) = 1/(2)^2;
        Q(11, 11) = 5000; % integral states
        Q(12, 12) = 5000; % integral states

        % R matrix
        R = zeros(nInputs);
        R(1, 1) = 1/(0.25)^2;
        R(2, 2) = 1/(0.5)^2;
        R(3, 3) = 1/(3)^2;
        R(4, 4) = 1/(10)^2;
        
        % measure first two joint angles
        C = zeros(2, nStates);
        C(1, 2) = 1;
        C(2, 3) = 1;

        A_aug = [A zeros(10, 2); -C zeros(2, 2)];
        B_aug = [B; zeros(2, 4)];

        K = lqr(A_aug, B_aug, Q, R); % lqr gains
        
%         des_poles = eig(A_aug-B_aug*K);
%         des_poles(7:8) = [ -0.9 + 0.8797i; -0.9- 0.8797i]; 
%         K = place(A_aug, B_aug,des_poles);
%         sysOL = ss(A,B,eye(nState),0);
%         sysOL_lqi = ss(A_aug,B_aug,eye(nState+2),0);
        sysCL_lqi = ss(A_aug-B_aug*K,B_aug,eye(nStates+2),0);
        
        figure;
        pzplot(sysCL_lqi);
        
        p_ctrl_fl.lqi_gain = K;
        kr_q = zeros(4);%-inv(C([1:4], 1:10)*inv(A - B*K)*B);
        p_ctrl_fl.kr = kr_q; %-inv(C*inv(A - B*K_fl)*B); % scale the reference appropriately

    %     Kr = -inv(C*inv(A - B*K)*B);

        p_ctrl_fl.zf = [zf; 0; 0]; % desired final state, augmented by two integral states

        control_law = @(t, z) control_law_feedback_linearization_with_ball_lqi(t, z, p_model, p_ctrl_fl);
        
        % ------------------ Integral dynamics -------------------%
        Cint = zeros(2, 12);
        Cint(1, 2) = 1;
        Cint(2, 3) = 1;
        
        int_dynamics  = @(z) (zf(2:3) - Cint*z);
        
        % ------------------ Observer Dynamics -----------------------%
        
%          obs_poles = eig(A_aug-B_aug*K)*5;
%         obs_poles = [-26.2458+26.6328i,...
%             -26.2458-26.6328i,-19.4203+26.4199i,...
%             -19.4203-26.4199i,-24.8976+0.0000i,-4.4256+4.3987i, -4.4256-4.3987i,...
%             -9.2173+14.3300i,-9.2173-14.3300i,-15]*2;
%          -10.7147 +24.0728i
     obs_poles = [-10.7147 + 24.0728i; -10.7147 - 24.0728i;
                  -23.7867 + 9.7110i; -23.7867 - 9.7110i;
                  -9.9250 + 17.3362i; -9.9250 - 17.3362i;
                  -5.6026 + 9.6764i; -5.6026 - 9.6764i;
                  -20.1056 + 0.0000i; -11.1640 + 0.0000i];

 %  obs_poles = obs_poles(1:10);
%         Cob_aug = [Cob zeros(5,2)];
        % pole placement (replace by Kalman later)
%         G=B(:,1);
%         G= eye(10,4);
%           G = zeros(10,1);G(1) =1;
%         N=zeros(4,5)
%         L = lqe(A,G,Cob,Wi,Wo);
%         sys = ss(A,[G],Cob,0);
%         L = kalman(sys,Wi(1,1),Wo)
%         [kest,L_opt,SIGMA] = kalman(ss(A,eye(10,4),Cob,[]),eye(4,4),Wo);
        
        L = place(A', C_ob', obs_poles)';
        p_obsv.L = L;
        p_obsv.A = A;
        p_obsv.B = B;
        p_obsv.C = C_ob;
        p_obsv.zf = zf;%p_ctrl_fl.zf; % augmented desired 
        observer_dynamics = @(y, x_hat, u) obsv_feedback_linearization_with_ball(y, x_hat, u, p_model, p_obsv);
                
    case 'operational_space_fb_lin'
    % design lqr for feedback-linearized system in operational space
        nDof = 2;
        nStates = nDof*2;
        
        A = [zeros(nDof) eye(nDof); zeros(nDof) zeros(nDof)];
        B = [zeros(nDof); eye(nDof)];
        C = eye(nStates); % full state feedback for identity matrix, only measure (x,y)
        assert(rank(ctrb(A, B)) == length(A), 'System is not controllable\n'); % check controllability

        % Q matrix
        Q = eye(nStates);
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

    case 'standard_lqr'
    % design standard lqr towards desired final state
        nDof = 5;
        nStates = nDof*2;
        nInputs = 4;
        
%         zf = zf(1:8); % ignore ball states

        % tangent linearization of dynamics about final desired position
        u_equi = Grav_arm(zf, p_model); % inputs at equilibrium = gravity - for non-zero set point
        
        [A_lin, B_lin] = linearize_dynamics(zf, u_equi, p_model);
        assert(rank(ctrb(A_lin, B_lin)) == length(A_lin), 'System is not controllable'); % check controllability
        assert(rank(obsv(A_lin, C_ob)) == length(A_lin), 'System is not observable'); % check observability
        fprintf(['Condition number of controllability  matrix: ', num2str(cond(ctrb(A_lin, B_lin))), '\n'])
        fprintf(['Condition number of observability matrix: ', num2str(cond(obsv(A_lin, C_ob))), '\n'])

        % Q matrix
        Q = zeros(nStates);
        Q(1, 1) = 1/(0.05)^2; % max cart position
        Q(2, 2) = 1/(deg2rad(10))^2; % max joint 1 angle
        Q(3, 3) = 1/(deg2rad(10))^2; % max joint 2 angle
        Q(4, 4) = 1/(deg2rad(3))^2; % max joint 3 angle
        Q(5, 5) = 0; % max cart velocity 
        Q(6, 6) = 1/(5)^2; % max joint 1 velocity
        Q(7, 7) = 1/(2.5)^2; % max joint 2 velocity
        Q(8, 8) = 1/(1)^2; % max joint 3 velocity
        Q(9, 9) = 1/(0.02)^2; % ball position deviation
        Q(10, 10) = 1/(3)^2; % ball velocity deviation

        % R matrix
        R = zeros(nInputs);
        R(1, 1) = 0.01; % cart force limit
        R(2, 2) = 1/(4)^2; % torque 1 limit
        R(3, 3) = 1/(2)^2; % torque 2 limit
        R(4, 4) = 1/(0.35)^2; % torque 3 limit

        K_lqr = lqr(A_lin, B_lin, Q, R);
        p_ctrl_lqr.lqr_gain = K_lqr;
        p_ctrl_lqr.zf = zf; % desired final state
        p_ctrl_lqr.kr = zeros(nInputs, nStates); % set reference input to 0
%         p_ctrl_lqr.kr = -(C_lin(1:4, :)/inv(A_lin - B_lin*K_lqr)*B_lin)*0;

        control_law = @(t, z) (control_law_standard_lqr(t, z, p_model, p_ctrl_lqr) + u_equi);
        
    %--------------------------- observer design ---------------------------- %
        cont_poles = eig(A_lin - B_lin*K_lqr);
        obs_poles = cont_poles*3; % arbitrary decision - poles 3x faster than system
     
        sys = ss(A_lin - B_lin*K_lqr, B_lin, C_ob, []);
        figure(15)
        pzmap(sys)
        grid on
        
        % pole placement (replace by Kalman later)
        L = place(A_lin', C_ob', obs_poles)';
        p_obsv.L = L;
        p_obsv.A = A_lin;
        p_obsv.B = B_lin;
        p_obsv.C = C_ob;
        p_obsv.zf = zf;
        p_obsv.u_equi = u_equi;
        observer_dynamics = @(y, x_hat, u) obsv_standard_lqr(y, x_hat, u, p_model, p_obsv);

    otherwise
        error('Choose controller from available options')
end