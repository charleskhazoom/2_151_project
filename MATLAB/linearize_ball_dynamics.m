function [A_lin, B_lin] = linearize_ball_dynamics(z_equi, u_equi, p)
% linearize_ball_dynamics: returns linearized system and input weighting
% matrices about given equilibrium state and inputs for ball dynamics only
%
% INPUTS
% z_equi: equilibrium state
% u_equi: equilibrium input
% p: parameters
%
% OUTPUTS
% A_lin: linearized system matrix
% B_lin: linearized input weighting matrix

    n_states = length(z_equi);
    n_inputs = length(u_equi);
    t = 0; % dummy assignation of t variable
    dz_equi = ball_dynamics(t, z_equi, u_equi, p); % equilibrium dz
    
    ep = 1e-6; % tolerance
    delStateMat = eye(n_states);
    delInputMat = eye(n_inputs);
    A_lin = zeros(2,n_states);
    B_lin = zeros(2, n_inputs);

    % linearize A matrix
    for i = 1:n_states
        z_dev = z_equi + ep*delStateMat(:, i); % small perturbation of state
        dz_dev = ball_dynamics(t, z_dev, u_equi, p); % rate of change of state from this position
        A_lin(:, i) = (dz_dev - dz_equi)/ep; % tangent to this perturbation
    end

    % linearize B matrix
    for i = 1:n_inputs
        u_dev = u_equi + ep*delInputMat(:, i); % small perturbation of input
        dz_dev = ball_dynamics(t, z_equi, u_dev, p); % rate of change of state from this input
        B_lin(:, i) = (dz_dev - dz_equi)/ep; % tanget to this perturbation
    end
end