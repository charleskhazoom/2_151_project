function b = b_arm(in1,in2,in3)
%B_ARM
%    B = B_ARM(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    24-Nov-2020 19:31:41

F_x = in2(1,:);
dth1 = in1(6,:);
dth2 = in1(7,:);
dx = in1(5,:);
g = in3(10,:);
l_1 = in3(7,:);
l_2 = in3(8,:);
m1 = in3(2,:);
m2 = in3(3,:);
m3 = in3(4,:);
tau1 = in2(2,:);
tau2 = in2(3,:);
tau3 = in2(4,:);
th1 = in1(2,:);
th2 = in1(3,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = l_1.^2;
t10 = -dx;
t6 = l_1.*t2;
t7 = cos(t4);
t8 = l_1.*t3;
t9 = sin(t4);
t11 = l_2.*t7;
t12 = l_2.*t9;
t17 = (dth1.*t8)./2.0;
t13 = dth1.*t11;
t14 = dth2.*t11;
t15 = dth1.*t12;
t16 = dth2.*t12;
t19 = -t17;
t20 = t11./2.0;
t21 = t12./2.0;
t26 = t6+t11;
t27 = t8+t12;
t18 = t14.*2.0;
t22 = t13./2.0;
t23 = t14./2.0;
t24 = t15./2.0;
t25 = t16./2.0;
t28 = dx+t19;
t29 = dth1.*t26;
t30 = dth1.*t27;
t31 = t6+t20;
t32 = t8+t21;
t33 = t13+t14;
t34 = t15+t16;
t35 = dth1.*t31;
t36 = dth1.*t32;
t37 = t22+t23;
t38 = t24+t25;
t39 = t14+t29;
t40 = t16+t30;
t41 = t10+t40;
t42 = t23+t35;
t43 = t25+t36;
t44 = t12.*t39.*2.0;
t45 = t10+t43;
t46 = t11.*t41.*2.0;
t47 = t12.*t42;
t48 = -t46;
t49 = t11.*t45;
t50 = -t49;
b = [F_x+dth1.*((m3.*(t18+t29.*2.0))./2.0+(m2.*(t14+t35.*2.0))./2.0+(dth1.*m1.*t6)./2.0)+dth2.*((m2.*t33)./2.0+(m3.*(t13.*2.0+t18))./2.0);tau1+dth1.*((m1.*(t6.*t28+(dth1.*t2.*t3.*t5)./2.0))./2.0+(m3.*(t26.*t40.*2.0-t26.*t41.*2.0))./2.0+(m2.*(t31.*t43.*2.0-t31.*t45.*2.0))./2.0)-(m1.*(dth1.*t6.*t28+(dth1.^2.*t2.*t3.*t5)./2.0))./2.0+dth2.*((m3.*(t44+t48+t26.*t34.*2.0-t27.*t33.*2.0))./2.0+(m2.*(t47+t50+t31.*t38.*2.0-t32.*t37.*2.0))./2.0)-(m3.*(t39.*t40.*2.0-t39.*t41.*2.0))./2.0-(m2.*(t42.*t43.*2.0-t42.*t45.*2.0))./2.0-(g.*m1.*t6)./2.0-g.*m3.*t26-g.*m2.*t31;tau2+dth2.*((m3.*(t44+t48+t11.*t34.*2.0-t12.*t33.*2.0))./2.0+(m2.*(t47+t50+t11.*t38-t12.*t37))./2.0)-(m3.*(t34.*t39.*2.0-t33.*t41.*2.0))./2.0-(m2.*(t38.*t42.*2.0-t37.*t45.*2.0))./2.0-dth1.*((m3.*(t46-t11.*t40.*2.0))./2.0+(m2.*(t49-t11.*t43))./2.0)-(g.*m2.*t11)./2.0-g.*m3.*t11;tau3];
