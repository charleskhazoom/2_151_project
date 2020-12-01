function u = control_law_feedback_linearization(t, z, p, p_ctrl)
% control_law_feedback_linearization: calculate next step control input
% using a feedback linearization scheme
%
% INPUTS
% t: time
% z: state
% p: parameters
% p_ctrl: controller parameters
%
% OUTPUT
% u: Control output

    z = z(1:8); % ignore ball states
    
    % unpack controller parameters
    K_lqr = p_ctrl.lqr_gain;
    kr = p_ctrl.kr; % scale the reference appropriately
    zf = p_ctrl.zf; % desired final state
    
    M = A_arm(z, p); % Mass matrix
    V = Corr_arm(z, p); % Coriolis matrix
    G = Grav_arm(z, p); % Gravity matrix
    
    r = zf(1:4); % desired new equilibrium point 
    w = kr*r - K_lqr*z(1:8); % do lqr here 
    u = V + G + M*w; % control input cancels nonlinear dynamics, adds state feedback
end