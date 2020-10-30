function Corr_Joint_Sp = Corr_arm(in1,in2)
%CORR_ARM
%    CORR_JOINT_SP = CORR_ARM(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    30-Oct-2020 16:09:10

dth1 = in1(6,:);
dth2 = in1(7,:);
dth3 = in1(8,:);
l_2 = in2(7,:);
l_3 = in2(8,:);
l_4 = in2(9,:);
m2 = in2(2,:);
m3 = in2(3,:);
m4 = in2(4,:);
th1 = in1(2,:);
th2 = in1(3,:);
th3 = in1(4,:);
t2 = cos(th1);
t3 = sin(th2);
t4 = sin(th3);
t5 = th1+th2;
t6 = th2+th3;
t7 = dth1.^2;
t8 = dth2.^2;
t9 = dth3.^2;
t10 = cos(t5);
t11 = sin(t6);
t12 = t5+th3;
t14 = dth1.*dth3.*l_3.*l_4.*m4.*t4;
t15 = dth2.*dth3.*l_3.*l_4.*m4.*t4;
t18 = (l_3.*l_4.*m4.*t4.*t9)./2.0;
t13 = cos(t12);
t16 = -t14;
t17 = -t15;
t19 = -t18;
Corr_Joint_Sp = [l_2.*m2.*t2.*t7.*(-1.0./2.0)-l_2.*m3.*t2.*t7-l_2.*m4.*t2.*t7-(l_3.*m3.*t7.*t10)./2.0-(l_3.*m3.*t8.*t10)./2.0-l_3.*m4.*t7.*t10-l_3.*m4.*t8.*t10-(l_4.*m4.*t7.*t13)./2.0-(l_4.*m4.*t8.*t13)./2.0-(l_4.*m4.*t9.*t13)./2.0-dth1.*dth2.*l_3.*m3.*t10-dth1.*dth2.*l_3.*m4.*t10.*2.0-dth1.*dth2.*l_4.*m4.*t13-dth1.*dth3.*l_4.*m4.*t13-dth2.*dth3.*l_4.*m4.*t13;t16+t17+t19-(l_2.*l_3.*m3.*t3.*t8)./2.0-l_2.*l_3.*m4.*t3.*t8-(l_2.*l_4.*m4.*t8.*t11)./2.0-(l_2.*l_4.*m4.*t9.*t11)./2.0-dth1.*dth2.*l_2.*l_3.*m3.*t3-dth1.*dth2.*l_2.*l_3.*m4.*t3.*2.0-dth1.*dth2.*l_2.*l_4.*m4.*t11-dth1.*dth3.*l_2.*l_4.*m4.*t11-dth2.*dth3.*l_2.*l_4.*m4.*t11;t16+t17+t19+(l_2.*l_3.*m3.*t3.*t7)./2.0+l_2.*l_3.*m4.*t3.*t7+(l_2.*l_4.*m4.*t7.*t11)./2.0;(l_4.*m4.*(l_3.*t4.*t7+l_3.*t4.*t8+l_2.*t7.*t11+dth1.*dth2.*l_3.*t4.*2.0))./2.0];
