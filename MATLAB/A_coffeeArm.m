function A = A_coffeeArm(in1,in2)
%A_COFFEEARM
%    A = A_COFFEEARM(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    30-Oct-2020 18:09:34

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
t5 = l_2.^2;
t6 = l_3.^2;
t7 = l_4.^2;
t8 = l_2.*t2;
t9 = cos(t4);
t10 = l_2.*t3;
t11 = sin(t4);
t12 = t4+th3;
t24 = (m3.*t6)./1.2e+1;
t25 = (m4.*t7)./1.2e+1;
t13 = cos(t12);
t14 = sin(t12);
t15 = t10.*2.0;
t16 = l_3.*t9;
t17 = l_3.*t11;
t20 = (m2.*t10)./2.0;
t18 = t17.*2.0;
t19 = l_4.*t14;
t21 = -t20;
t22 = t16./2.0;
t23 = t17./2.0;
t26 = (l_4.*t13)./2.0;
t30 = m3.*t17.*(-1.0./2.0);
t32 = t15+t17;
t27 = m3.*t23;
t28 = t19./2.0;
t31 = m4.*t19.*(-1.0./2.0);
t33 = t8+t22;
t34 = t10+t23;
t35 = t18+t19;
t36 = (m3.*t32)./2.0;
t37 = t16+t26;
t29 = m4.*t28;
t38 = t17+t28;
t39 = -t36;
t40 = t17.*t34;
t41 = t15+t35;
t42 = (m4.*t35)./2.0;
t43 = t16.*t33;
t44 = t8+t37;
t47 = l_4.*t13.*t37;
t45 = t10+t38;
t46 = -t42;
t48 = t19.*t38;
t49 = (m4.*t41)./2.0;
t51 = l_4.*t13.*t44;
t54 = t37.*t44.*2.0;
t56 = t40+t43;
t50 = -t49;
t52 = t19.*t45;
t53 = t30+t46;
t55 = t38.*t45.*2.0;
t57 = (m3.*t56)./2.0;
t58 = t47+t48;
t59 = t21+t39+t50;
t60 = (m4.*t58)./2.0;
t61 = t51+t52;
t65 = t54+t55;
t62 = (m4.*t61)./2.0;
t63 = t25+t60;
t66 = (m4.*t65)./2.0;
t64 = t25+t62;
t67 = t24+t25+t57+t66;
A = reshape([m1+m2+m3+m4,t59,t53,t31,t59,t24+t25+(m2.*t5)./1.2e+1+(m2.*((t2.^2.*t5)./2.0+(t3.^2.*t5)./2.0))./2.0+(m3.*(t33.^2.*2.0+t34.^2.*2.0))./2.0+(m4.*(t44.^2.*2.0+t45.^2.*2.0))./2.0,t67,t64,t53,t67,t24+t25+(m3.*((t6.*t9.^2)./2.0+(t6.*t11.^2)./2.0))./2.0+(m4.*(t37.^2.*2.0+t38.^2.*2.0))./2.0,t63,t31,t64,t63,t25+(m4.*((t7.*t13.^2)./2.0+(t7.*t14.^2)./2.0))./2.0],[4,4]);
