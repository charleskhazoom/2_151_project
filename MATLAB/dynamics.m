function dz = dynamics(t,z,u,p)
% calculate rate of change of state dz at time t from state z, input u, and
% system parameters p

    % Get mass matrix
    A = A_arm(z,p);
     
    % Get b = Q - V(q,qd) - G(q)
    b = b_arm(z,u,p);
    
    % Solve for qdd.
    qdd = A\b;
    dz = 0*z;
    
    % Form dz
    dz(1:4) = z(5:8);
    dz(5:8) = qdd;
end