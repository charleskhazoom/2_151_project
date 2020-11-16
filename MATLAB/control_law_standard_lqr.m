function u = control_law_standard_lqr(t, z, p, p_ctrl)
% find next step control input using standard lqr

    % extract control parameters
    K_lqr = p_ctrl.lqr_gain;
    K_r = p_ctrl.kr;
    zf = p_ctrl.zf;
    
    % don't want dz componenets of final state - just want positions
%     zf = zf(1:4);
    
    % u = r - K*x
    u = K_r*zf(1:4) - K_lqr*(z-zf);
    
end