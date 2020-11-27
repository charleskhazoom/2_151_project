function dz_ball = ball_dynamics(t,z,qdd,p)
        qdot = z(5:8);
        %% adds the ball state to the dynamics
        dz_ball(1) = z(10); 
        
        theta = z(2) + z(3) - 90/180*pi + z(4); % angle of plate
        
        jdot = jacobian_dot_endEffector(z, p);
        j = jacobian_endEffector(z, p);
        %xdd = jdot*dz(1:4) + j*qdd; % get linear acceleration using jacobians 
        xdd = jdot*qdot + j*qdd; % get linear acceleration using jacobians 
        
        a_x_plate = xdd(1)*cos(-theta) - xdd(2)*sin(-theta); % rotate into plate frame        
        dz_ball(2) = (1/(1 + 2/5))*(p(10)*sin(-theta) - a_x_plate - 1*z(10)); % p(10) = gravity and a damping term -k*v
        dz_ball = dz_ball';
end