function E = energy_coffeeArm(in1,in2)
%ENERGY_COFFEEARM
%    E = ENERGY_COFFEEARM(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    04-Nov-2020 23:12:31

dth1 = in1(6,:);
dth2 = in1(7,:);
dth3 = in1(8,:);
dx = in1(5,:);
g = in2(10,:);
l_2 = in2(7,:);
l_3 = in2(8,:);
l_4 = in2(9,:);
m1 = in2(1,:);
m2 = in2(2,:);
m3 = in2(3,:);
m4 = in2(4,:);
th1 = in1(2,:);
th2 = in1(3,:);
th3 = in1(4,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = dth1.^2;
t6 = l_2.^2;
t12 = -dx;
t7 = l_2.*t2;
t8 = cos(t4);
t9 = l_2.*t3;
t10 = sin(t4);
t11 = t4+th3;
t13 = cos(t11);
t14 = sin(t11);
t15 = l_3.*t8;
t16 = l_3.*t10;
t17 = t16./2.0;
t18 = (l_4.*t13)./2.0;
t19 = (l_4.*t14)./2.0;
t20 = t9+t17;
t21 = t9+t16+t19;
E = (m3.*(((dth2.*t15)./2.0+dth1.*(t7+t15./2.0)).^2+(t12+dth2.*t17+dth1.*t20).^2))./2.0+(m4.*((dth2.*(t15+t18)+dth3.*t18+dth1.*(t7+t15+t18)).^2+(t12+dth2.*(t16+t19)+dth1.*t21+dth3.*t19).^2))./2.0+(dx.^2.*m1)./2.0+(m2.*((dx-(dth1.*t9)./2.0).^2+(t2.^2.*t5.*t6)./4.0))./2.0+(l_3.^2.*m3.*(dth1+dth2).^2)./2.4e+1+(l_4.^2.*m4.*(dth1+dth2+dth3).^2)./2.4e+1+(g.*m2.*t9)./2.0+g.*m3.*t20+g.*m4.*t21+(m2.*t5.*t6)./2.4e+1;
