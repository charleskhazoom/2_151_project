function r4 = position_endEffector(in1,in2)
%POSITION_ENDEFFECTOR
%    R4 = POSITION_ENDEFFECTOR(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    30-Oct-2020 16:09:10

l_2 = in2(7,:);
l_3 = in2(8,:);
l_4 = in2(9,:);
th1 = in1(2,:);
th2 = in1(3,:);
th3 = in1(4,:);
x = in1(1,:);
t2 = th1+th2;
t3 = t2+th3;
r4 = [x+l_3.*cos(t2)+l_4.*cos(t3)+l_2.*cos(th1);l_3.*sin(t2)+l_4.*sin(t3)+l_2.*sin(th1)];