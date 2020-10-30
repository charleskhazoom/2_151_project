function dJ = jacobian_dot_endEffector(in1,in2)
%JACOBIAN_DOT_ENDEFFECTOR
%    DJ = JACOBIAN_DOT_ENDEFFECTOR(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    30-Oct-2020 16:09:10

dth1 = in1(6,:);
dth2 = in1(7,:);
dth3 = in1(8,:);
l_2 = in2(7,:);
l_3 = in2(8,:);
l_4 = in2(9,:);
th1 = in1(2,:);
th2 = in1(3,:);
th3 = in1(4,:);
t2 = th1+th2;
t3 = cos(t2);
t4 = sin(t2);
t5 = t2+th3;
t6 = cos(t5);
t7 = sin(t5);
t8 = l_3.*t3;
t9 = l_3.*t4;
t10 = l_4.*t6;
t11 = l_4.*t7;
dJ = reshape([0.0,0.0,-dth2.*(t8+t10)-dth3.*t10-dth1.*(t8+t10+l_2.*cos(th1)),-dth2.*(t9+t11)-dth3.*t11-dth1.*(t9+t11+l_2.*sin(th1))],[2,2]);
