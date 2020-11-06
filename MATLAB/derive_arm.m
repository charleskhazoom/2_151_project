clear all; clc

name = 'arm';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms t x th1 th2 th3 dx dth1 dth2 dth3 ddx ddth1 ddth2 ddth3 real
syms m_cart m1 m2 m3 h_cart l_cart l_1 l_2 l_3 g real % distance to CoM for links is half l 
syms F_x tau1 tau2 tau3 real
%syms Ir N real % can add in rotor inertia and gearing later

% Group them
q   = [x; th1  ; th2 ; th3];      % generalized coordinates
dq  = [dx; dth1 ; dth2; dth3];   % first time derivatives
ddq = [ddx; ddth1; ddth2; ddth3];  % second time derivatives
u   = [F_x; tau1 ; tau2; tau3];     % controls

p   = [m_cart m2 m3 m4 h_cart l_cart l_1 l_2 l_3 g]';        % parameters

I1 = m1*l_1^2/12;
I2 = m2*l_2^2/12;
I3 = m3*l_3^2/12;

% Generate Vectors and Derivativess
ihat = [1; 0; 0];
jhat = [0; 1; 0];

xhat = [1; 0; 0];
yhat = [0; 1; 0];
khat = cross(ihat,jhat);

e1hat =  cos(th1)*ihat + sin(th1)*jhat;
e2hat =  cos(th1+th2)*ihat + sin(th1+th2)*jhat;
e3hat =  cos(th1+th2+th3-pi/2)*ihat + sin(th1+th2+th3-pi/2)*jhat;

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

rCart = x*ihat;
r1 = rCart + l_1 * e1hat;
r2 = r1  + l_2 * e2hat;
r3_R = r2 + l_3/2 * e3hat;
r3_L = r2 - l_3/2 * e3hat;

r_mCart = rCart;
r_m1 = rCart + (l_1/2)*e1hat;
r_m2 = r1 + (l_2/2)*e2hat;
r_m3 = r2;

drCart = ddt(rCart);
dr1 = ddt(r1);
dr2 = ddt(r2);
dr3 = ddt(r3);

dr_m1 = ddt(r_m1);
dr_m2 = ddt(r_m2);
dr_m3 = ddt(r_m3);

% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces

omega1 = dth1;
omega2 = omega1 + dth2;
omega3 = omega2 + dth3;

T_cart = (1/2)*m_cart * dot(drCart,drCart);
T1 = (1/2)*m1 * dot(dr_m1,dr_m1) + (1/2) * I1 * omega1^2;
T2 = (1/2)*m2 * dot(dr_m2,dr_m2) + (1/2) * I2 * omega2^2;
T3 = (1/2)*m3 * dot(dr_m3,dr_m3) + (1/2) * I3 * omega3^2;

% T1r = (1/2)*Ir*(N*dth1)^2;
% T2r = (1/2)*Ir*(dth1 + N*dth2)^2; % rotor inertia energies

Vg1 = m1*g*r_m1(2);
Vg2 = m2*g*r_m2(2);
Vg3 = m3*g*r_m3(2);

T = simplify(T_cart + T1 + T2 + T3);
Vg = Vg1 + Vg2 + Vg3;

Q_f = F2Q(F_x*ihat,rCart);
Q_tau1 = M2Q(tau1*khat,dth1*khat);
Q_tau2 = M2Q(tau2*khat,dth2*khat); 
Q_tau3 = M2Q(tau3*khat,dth3*khat);

Q_tau = Q_tau1 + Q_tau2 + Q_tau3;

Q = Q_tau + Q_f;

% Assemble the array of cartesian coordinates of the key points
keypoints = [rCart(1:2) r1(1:2) r2(1:2) r3_R(1:2) r3_L(1:2)];

%% All the work is done!  Just turn the crank...
% Derive Energy Function and Equations of Motion
E = T+Vg;
L = T-Vg;
eom = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;

% Rearrange Equations of Motion
A = jacobian(eom,ddq);
b = A*ddq - eom;

% Equations of motion are
% eom = A *ddq + (coriolis term) + (gravitational term) - Q = 0
Mass_Joint_Sp = A;
Grav_Joint_Sp = simplify(jacobian(Vg, q)');
Corr_Joint_Sp = simplify( eom + Q - Grav_Joint_Sp - A*ddq);

% Compute foot jacobian
J = jacobian(r3,q);

% Compute ddt( J )
dJ= reshape( ddt(J(:)) , size(J) );

% Write Energy Function and Equations of Motion
z  = [q ; dq];
r3 = r3(1:2);
dr3= dr3(1:2);
J  = J(1:2,1:4);
dJ = dJ(1:2,1:4);
%%
matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z u p});
matlabFunction(E,'file',['energy_' name],'vars',{z p});
matlabFunction(r3,'file',['position_endEffector'],'vars',{z p});
matlabFunction(dr3,'file',['velocity_endEffector'],'vars',{z p});
matlabFunction(J ,'file',['jacobian_endEffector'],'vars',{z p});
matlabFunction(dJ ,'file',['jacobian_dot_endEffector'],'vars',{z p});

matlabFunction(Grav_Joint_Sp ,'file', ['Grav_arm'] ,'vars',{z p});
matlabFunction(Corr_Joint_Sp ,'file', ['Corr_arm']     ,'vars',{z p});
matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});




