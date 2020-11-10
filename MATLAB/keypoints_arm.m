function keypoints = keypoints_arm(in1,in2)
%KEYPOINTS_ARM
%    KEYPOINTS = KEYPOINTS_ARM(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    09-Nov-2020 22:55:31

l_1 = in2(7,:);
l_2 = in2(8,:);
l_3 = in2(9,:);
th1 = in1(2,:);
th2 = in1(3,:);
th3 = in1(4,:);
x = in1(1,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t9 = pi./2.0;
t5 = l_1.*t2;
t6 = cos(t4);
t7 = l_1.*t3;
t8 = sin(t4);
t12 = -t9;
t10 = l_2.*t6;
t11 = l_2.*t8;
t13 = t4+t12+th3;
t14 = cos(t13);
t15 = sin(t13);
t16 = (l_3.*t14)./2.0;
t17 = (l_3.*t15)./2.0;
keypoints = reshape([x,0.0,t5+x,t7,t5+t10+x,t7+t11,t5+t10+t16+x,t7+t11+t17,t5+t10-t16+x,t7+t11-t17],[2,5]);
