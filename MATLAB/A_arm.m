function A = A_arm(in1,in2)
%A_ARM
%    A = A_ARM(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    24-Nov-2020 19:31:32

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
t5 = l_1.^2;
t6 = l_2.^2;
t7 = l_3.^2;
t8 = l_1.*t2;
t9 = cos(t4);
t10 = l_1.*t3;
t11 = sin(t4);
t23 = (m2.*t6)./1.2e+1;
t24 = (m3.*t7)./1.2e+1;
t12 = t10.*2.0;
t13 = t9.^2;
t14 = t11.^2;
t15 = l_2.*t9;
t16 = l_2.*t11;
t19 = (m1.*t10)./2.0;
t17 = t16.*2.0;
t18 = m3.*t16;
t20 = -t19;
t21 = t15./2.0;
t22 = t16./2.0;
t27 = t8+t15;
t28 = t10+t16;
t29 = m2.*t16.*(-1.0./2.0);
t30 = t12+t16;
t25 = -t18;
t26 = m2.*t22;
t31 = t12+t17;
t32 = t8+t21;
t33 = t10+t22;
t34 = (m2.*t30)./2.0;
t35 = t17.*t28;
t38 = t15.*t27.*2.0;
t36 = -t34;
t37 = (m3.*t31)./2.0;
t40 = t16.*t33;
t41 = t15.*t32;
t42 = t25+t29;
t43 = t35+t38;
t39 = -t37;
t44 = t40+t41;
t45 = (m3.*t43)./2.0;
t46 = (m2.*t44)./2.0;
t47 = t20+t36+t39;
t48 = t23+t24+t45+t46;
A = reshape([m1+m2+m3+m_cart,t47,t42,0.0,t47,t23+t24+(m1.*t5)./1.2e+1+(m1.*((t2.^2.*t5)./2.0+(t3.^2.*t5)./2.0))./2.0+(m3.*(t27.^2.*2.0+t28.^2.*2.0))./2.0+(m2.*(t32.^2.*2.0+t33.^2.*2.0))./2.0,t48,t24,t42,t48,t23+t24+(m3.*(t6.*t13.*2.0+t6.*t14.*2.0))./2.0+(m2.*((t6.*t13)./2.0+(t6.*t14)./2.0))./2.0,t24,0.0,t24,t24,t24],[4,4]);
