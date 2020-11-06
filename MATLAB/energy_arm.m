function E = energy_arm(in1,in2)
%ENERGY_ARM
%    E = ENERGY_ARM(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    05-Nov-2020 23:55:47

dth1 = in1(6,:);
dth2 = in1(7,:);
dth3 = in1(8,:);
dx = in1(5,:);
g = in2(10,:);
l_1 = in2(7,:);
l_2 = in2(8,:);
l_3 = in2(9,:);
m1 = in2(2,:);
m2 = in2(3,:);
m3 = in2(4,:);
m_cart = in2(1,:);
th1 = in1(2,:);
th2 = in1(3,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = dth1.^2;
t6 = l_1.^2;
t11 = -dx;
t7 = l_1.*t2;
t8 = cos(t4);
t9 = l_1.*t3;
t10 = sin(t4);
t12 = l_2.*t10;
t13 = t12./2.0;
t14 = t9+t12;
t15 = t9+t13;
E = (m3.*((dth1.*(t7+l_2.*t8)+dth2.*l_2.*t8).^2+(t11+dth2.*t12+dth1.*t14).^2))./2.0+(m2.*((dth1.*(t7+(l_2.*t8)./2.0)+(dth2.*l_2.*t8)./2.0).^2+(t11+dth2.*t13+dth1.*t15).^2))./2.0+(dx.^2.*m_cart)./2.0+(m1.*((dx-(dth1.*t9)./2.0).^2+(t2.^2.*t5.*t6)./4.0))./2.0+(l_2.^2.*m2.*(dth1+dth2).^2)./2.4e+1+(l_3.^2.*m3.*(dth1+dth2+dth3).^2)./2.4e+1+(g.*m1.*t9)./2.0+g.*m2.*t15+g.*m3.*t14+(m1.*t5.*t6)./2.4e+1;
