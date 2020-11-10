function config = invKin_arm(pose,p,guess)
    % INPUTS:
    % pose: [2 by 1] task space position of end effector
    % p: [10 by 1] parameters
    % guess: [4 by 1] initial guess for ik solver
    
    % OUTPUTS:
    % config: [4 by 1] joint configuration for pose
    
    % Define variables for time, generalized coordinates + derivatives, controls, and parameters 
    syms t th1 th2 real
    q = [th1 th2]';
    
    % unpack parameters
    l_1 = p(7);
    l_2 = p(8);
    l_3 = p(9);
    
    x = pose(1);
    th3 = 0;
    
    % Generate Vectors and Derivativess
    ihat = [1; 0; 0];
    jhat = [0; 1; 0];
    khat = cross(ihat,jhat);

    e1hat =  cos(th1)*ihat + sin(th1)*jhat;
    e2hat =  cos(th1+th2)*ihat + sin(th1+th2)*jhat;
    e3hat =  cos(th1+th2+th3-pi/2)*ihat + sin(th1+th2+th3-pi/2)*jhat;
    
    rCart = x*ihat;
    r1 = rCart + l_1 * e1hat;
    r2 = r1  + l_2 * e2hat;
    r3 = r2;
    
    eqn_x = r3(1) == pose(1);
    eqn_y = r3(2) == pose(2);
    
    soln = vpasolve([eqn_x eqn_y],[th1 th2],guess);
    th3 = pi/2-soln.th1-soln.th2;
    
    config = [x soln.th1 soln.th2 th3]';
end