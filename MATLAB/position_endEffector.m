function r3 = position_endEffector(in1,in2)
%POSITION_ENDEFFECTOR
%    R3 = POSITION_ENDEFFECTOR(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    24-Nov-2020 19:31:41

l_1 = in2(7,:);
l_2 = in2(8,:);
th1 = in1(2,:);
th2 = in1(3,:);
x = in1(1,:);
t2 = th1+th2;
r3 = [x+l_2.*cos(t2)+l_1.*cos(th1);l_2.*sin(t2)+l_1.*sin(th1)];
