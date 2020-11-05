function b = b_coffeeArm(in1,in2,in3)
%B_COFFEEARM
%    B = B_COFFEEARM(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    04-Nov-2020 23:12:30

F_x = in2(1,:);
dth1 = in1(6,:);
dth2 = in1(7,:);
dth3 = in1(8,:);
dx = in1(5,:);
g = in3(10,:);
l_2 = in3(7,:);
l_3 = in3(8,:);
l_4 = in3(9,:);
m2 = in3(2,:);
m3 = in3(3,:);
m4 = in3(4,:);
tau1 = in2(2,:);
tau2 = in2(3,:);
tau3 = in2(4,:);
th1 = in1(2,:);
th2 = in1(3,:);
th3 = in1(4,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = l_2.^2;
t11 = -dx;
t6 = l_2.*t2;
t7 = cos(t4);
t8 = l_2.*t3;
t9 = sin(t4);
t10 = t4+th3;
t12 = cos(t10);
t13 = sin(t10);
t14 = l_3.*t7;
t15 = l_3.*t9;
t17 = (dth1.*t8)./2.0;
t16 = dth2.*t14;
t18 = dth3.*l_4.*t12;
t19 = -t17;
t20 = t14./2.0;
t21 = t15./2.0;
t24 = (l_4.*t12)./2.0;
t27 = (l_4.*t13)./2.0;
t22 = dth1.*t20;
t23 = t16./2.0;
t25 = dth1.*t21;
t26 = dth2.*t21;
t28 = dth1.*t24;
t29 = dth2.*t24;
t30 = t18./2.0;
t31 = dth1.*t27;
t32 = dth2.*t27;
t33 = dth3.*t27;
t34 = dx+t19;
t35 = t6+t20;
t36 = t8+t21;
t39 = t14+t24;
t40 = t15+t27;
t37 = dth1.*t35;
t38 = dth1.*t36;
t41 = dth1.*t39;
t42 = dth2.*t39;
t43 = dth1.*t40;
t44 = dth2.*t40;
t46 = t22+t23;
t47 = t25+t26;
t48 = t6+t39;
t49 = t8+t40;
t56 = t28+t29+t30;
t57 = t31+t32+t33;
t45 = t42.*2.0;
t50 = dth1.*t48;
t51 = dth1.*t49;
t52 = t23+t37;
t53 = t26+t38;
t60 = t33+t43+t44;
t61 = t30+t41+t42;
t54 = t11+t53;
t55 = t15.*t52;
t62 = t30+t42+t50;
t63 = t33+t44+t51;
t58 = t14.*t54;
t64 = t11+t63;
t65 = l_4.*t13.*t62;
t68 = t40.*t62.*2.0;
t59 = -t58;
t66 = -t65;
t67 = l_4.*t12.*t64;
t69 = t39.*t64.*2.0;
t70 = -t69;
b = [F_x+dth1.*((m3.*(t16+t37.*2.0))./2.0+(m4.*(t18+t45+t50.*2.0))./2.0+(dth1.*m2.*t6)./2.0)+dth2.*((m3.*(t16+dth1.*t14))./2.0+(m4.*(t18+t41.*2.0+t45))./2.0)+(dth3.*m4.*(t18+dth1.*l_4.*t12+dth2.*l_4.*t12))./2.0;tau1+dth1.*((m2.*(t6.*t34+(dth1.*t2.*t3.*t5)./2.0))./2.0+(m3.*(t35.*t53.*2.0-t35.*t54.*2.0))./2.0+(m4.*(t48.*t63.*2.0-t48.*t64.*2.0))./2.0)-(m2.*(dth1.*t6.*t34+(dth1.^2.*t2.*t3.*t5)./2.0))./2.0+dth2.*((m3.*(t55+t59+t35.*t47.*2.0-t36.*t46.*2.0))./2.0+(m4.*(t68+t70+t48.*t60.*2.0-t49.*t61.*2.0))./2.0)-(m3.*(t52.*t53.*2.0-t52.*t54.*2.0))./2.0-(m4.*(t62.*t63.*2.0-t62.*t64.*2.0))./2.0+(dth3.*m4.*(t65-t67+t48.*t57.*2.0-t49.*t56.*2.0))./2.0-(g.*m2.*t6)./2.0-g.*m3.*t35-g.*m4.*t48;tau2+dth2.*((m3.*(t55+t59+t14.*t47-t15.*t46))./2.0+(m4.*(t68+t70+t39.*t60.*2.0-t40.*t61.*2.0))./2.0)-(m3.*(t47.*t52.*2.0-t46.*t54.*2.0))./2.0-(m4.*(t60.*t62.*2.0-t61.*t64.*2.0))./2.0-dth1.*((m3.*(t58-t14.*t53))./2.0+(m4.*(t69-t39.*t63.*2.0))./2.0)+(dth3.*m4.*(t65-t67+t39.*t57.*2.0-t40.*t56.*2.0))./2.0-(g.*m3.*t14)./2.0-g.*m4.*t39;tau3-(m4.*(t57.*t62.*2.0-t56.*t64.*2.0))./2.0+(dth3.*m4.*(t65-t67+l_4.*t12.*t57-l_4.*t13.*t56))./2.0+(dth2.*m4.*(t65-t67+l_4.*t12.*t60-l_4.*t13.*t61))./2.0-(dth1.*m4.*(t67-l_4.*t12.*t63))./2.0-(g.*l_4.*m4.*t12)./2.0];
