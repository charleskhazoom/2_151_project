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