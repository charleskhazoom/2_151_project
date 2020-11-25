function u = control_law_standard_lqr(t, z, p, p_ctrl)
% control_law_standard_lqr: calculate next step control input using
% standard LQR control scheme
%
% INPUTS
% t: time
% z: state
% p: parameters
% p_ctrl: controller parameters
%
% OUTPUT
% u: Control output

    % extract controller parameters
    K_lqr = p_ctrl.lqr_gain;
    K_r = p_ctrl.kr;
    zf = p_ctrl.zf;
    
    % u = r - K*x
    u = K_r*zf(1:4) - K_lqr*(z(1:8)-zf(1:8));
end