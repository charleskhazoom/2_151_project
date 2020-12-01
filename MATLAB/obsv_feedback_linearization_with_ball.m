function dz_hat = obsv_feedback_linearization_with_ball(y, x_hat, u, p, p_obsv)
% obsv_feedback_linearization_with_ball: calculate observer dynamics for the
% feedback linearization controller with ball states
%
% INPUTS
% y: measured output
% x_hat: estimated state
% u: control input
% p: parameters
% p_obsv: observer parameters
%
% OUTPUTS
% dz_hat: rate of change of observer estimate

    A = p_obsv.A;
    B = p_obsv.B;
    C = p_obsv.C;
    L = p_obsv.L;
    zf = p_obsv.zf;

    % undo feedback linearization
    M = A_arm(x_hat, p); % Mass matrix
    V = Corr_arm(x_hat, p); % Coriolis matrix
    G = Grav_arm(x_hat, p); % Gravity matrix
    w = M\(u - V - G);

    % compute observer dynamics to be integrated
    dz_hat = (A - L*C)*(x_hat - zf) + B*w + L*(y - C*zf);
end