function u = control_law_fb_lin_op_sp(t, z, p, p_ctrl)
% control_law_fb_lin_op_sp: calculate next step control input according to
% feedback linearization scheme with operational state
%
% INPUTS
% t: time
% z: state
% p: parameters
% p_ctrl: controller parameters
%
% OUTPUT
% u: Control output

    K_lqr = p_ctrl.lqr_gain;
    kr = p_ctrl.kr; % scale the reference appropriately
    r = p_ctrl.r; % desired final state

    J = jacobian_endEffector(z, p);
    dJ = jacobian_dot_endEffector(z, p);
    
    M = A_arm(z, p); % system matrix
    V = Corr_arm(z, p); % Coriolis matrix
    G = Grav_arm(z, p); % Gravitational terms
    
    inv_M = inv(M);
    lambda = inv(J*inv_M*J');
    
    mu = lambda*J*inv_M*V - lambda*dJ*z(5:8);
    rho = lambda*J*inv_M*G;
    
    rE = position_endEffector(z, p); % position of end effector
    vE = velocity_endEffector(z, p); % velocity of end effector

    % end effector state
    operational_state = [rE; vE]; % just [x; y; dx; dy]

    w = kr*r - K_lqr*operational_state; % do lqr here 
    
    % damping_matrix = eye(4)*0.01;
    u = J'*(mu + rho + lambda\w);% - damping_matrix*z(5:8);
    
end