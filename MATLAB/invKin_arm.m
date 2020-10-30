function config = invKin_arm(pose,p,guess)
    % INPUTS:
    % pose: [2 by 1] task space position of end effector
    % p: [10 by 1] parameters
    % guess: [4 by 1] initial guess for ik solver
    
    % OUTPUTS:
    % config: [4 by 1] joint configuration for pose
    
    % Define variables for time, generalized coordinates + derivatives, controls, and parameters 
    syms t x th1 th2 th3 real
    q = [x th1 th2 th3]';
    
    % unpack parameters
    l_2 = p(7);
    l_3 = p(8);
    l_4 = p(9);
    
    % Generate Vectors and Derivativess
    ihat = [1; 0; 0];
    jhat = [0; 1; 0];
    khat = cross(ihat,jhat);

    e2hat =  cos(th1)*ihat + sin(th1)*jhat;
    e3hat =  cos(th1+th2)*ihat + sin(th1+th2)*jhat;
    e4hat =  cos(th1+th2+th3)*ihat + sin(th1+th2+th3)*jhat;
    
    rCart = x*ihat;
    r2 = rCart + l_2 * e2hat;
    r3 = r2  + l_3 * e3hat;
    r4 = r3  + l_4 * e4hat;
    
    eqn_x = r4(1) == pose(1);
    eqn_y = r4(2) == pose(2);
    soln = vpasolve([eqn_x eqn_y],[x th1 th2 th3],guess);
    config = [soln.x soln.th1 soln.th2 soln.th3]';
end