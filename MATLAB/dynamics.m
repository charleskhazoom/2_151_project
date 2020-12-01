function dz = dynamics(t, z, u, p)
% dynamics: calculate rate of change of state dz at time t from state z, input u, and
% system parameters p

% INPUTS
% t: time
% z: state
% u: control input
% p: parameters
%
% OUTPUTS
% dz: rate of change of state (system dynamics)
    
    % Get mass matrix
    A = A_arm(z, p);
     
    % Get b = Q - V(q,qd) - G(q)
    % non-acceleration terms
    b = b_arm(z, u, p);
    
    % Solve for qdd (acceleration term)
    qdd = A\b;
    dz = 0*z; % sets size of rate of change of states
    
    % Form dz
    dz(1:4) = z(5:8);
    dz(5:8) = qdd;
    
    % for the ball
    if (length(z) == 8)
    elseif (length(z) == 10)
        dz_ball = ball_dynamics(t, z, qdd, p);
        dz(9) = dz_ball(1);
        dz(10) = dz_ball(2);
    end
    
end