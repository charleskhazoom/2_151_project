function dz_hat = obsv_standard_lqr(y, x_hat, u, p, p_obsv)
% obsv_standard_lqr: calculate observer dynamics for the standard lqr
% controller
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

    % compute observer dynamics to be integrated
    dz_hat = (A - L*C)*x_hat + B*u + L*y;
end