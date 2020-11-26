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
    
    % Get system matrix
    A = A_arm(z,p);
     
    % Get b = Q - V(q,qd) - G(q)
    % Input weighting matrix
    b = b_arm(z,u,p);
    
    % Solve for qdd (acceleration term)
    qdd = A\b;
    dz = 0*z; % sets size of rate of change of states
    
    % Form dz
    dz(1:4) = z(5:8);
    dz(5:8) = qdd;
    
    % for the ball
    if (length(z) == 8)
    elseif (length(z) == 10)
        dz(9) = z(10); 
        
        theta = z(2) + z(3) - 90/180*pi + z(4); % angle of plate
        
        jdot = jacobian_dot_endEffector(z, p);
        j = jacobian_endEffector(z, p);
        xdd = jdot*dz(1:4) + j*qdd; % get linear acceleration using jacobians 
        
        a_x_plate = xdd(1)*cos(-theta) - xdd(2)*sin(-theta); % rotate into plate frame        
        dz(10) = (1/(1 + 2/5))*(p(10)*sin(-theta) - a_x_plate - 10*z(10)); % p(10) = gravity and a damping term -k*v
    end
end