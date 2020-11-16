function [A_lin, B_lin] =  linearize_dynamics(z_equi,u_equi,p)
% linearize the system dynamics about equilibrium state z_equi, input
% u_equi, and system parameters p

    n_states = length(z_equi);
    n_inputs = length(u_equi);
    t =0; % dummy assignation of t variable
    dz_equi = dynamics(t,z_equi,u_equi,p);
%     zout0 = zout(:,end);
    
    %statef0
    ep = 1e-6;
    delStateMat = eye(n_states);
    A_lin = zeros(n_states,n_states);
    B_lin = zeros(n_states,n_inputs);

    for i=1:n_states
        z_dev = z_equi + ep*delStateMat(:,i);
        dz_dev = dynamics(t,z_dev,u_equi,p);
        A_lin(:,i) = (dz_dev - dz_equi)/ep;
    end


    del_u_mat= eye(n_inputs);
    for i=1:n_inputs
        u_dev = u_equi+ep*del_u_mat(:,i);
        dz_dev = dynamics(t,z_equi,u_dev,p);
        B_lin(:,i) = (dz_dev - dz_equi)/ep;
    end
end