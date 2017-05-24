%% LQR
C2 = eye(6);
D2 = zeros(6,3);
sys2 = ss(A,B,C2,D2);
sysdeul2 = c2d(sys, Ts, 'forward');
Ad2 = sysdeul2.A;
Bd2 = sysdeul2.B;
Cd2 = sysdeul2.C;
Dd2 = sysdeul2.D;
Q1= eye(6);
R1 = eye(3) * 10^-4;
R =  eye(3) * 10^-4;
K_lqr1 = lqr(sysdrect2,Q1,R1,0);
Q2 = eye(6);
Q2(4:6,4:6) = 0;
K_lqr2 = lqr(sysdrect2,Q2,R1,0);
Q = eye(6);
Q(3,3) = 100;
K_lqr3 = lqr(Ad2,Bd2,Q,R,0);
x0 = [0;0;0;0;0;0];

uref = uopt;
xref = Y(:,2:end);
d1 = digits;
digits(3);
latexk1 = latex(sym(K_lqr1,'d'));
latexk2 = latex(sym(K_lqr2,'d'));
digits(d1)

%% Sim with LQR
close all;
Fmax = 200;
Mmax = 10;
Y2 = [xref(:,1)];
U2 = [];
for i = 1:N
    ufromx = K_lqr1*(Y2(:,i)-xref(:,i));
    u = uref(:,i) - ufromx;
    u = clipu(u, Fmax, Mmax);
    U2 = [U2,u];
    xnew = Ad2*Y2(:,i)+Bd2*u;
    ynew = xnew;
    Y2 = [Y2,ynew];
end

Y3 = [xref(:,1)];
U3 = [];
for i = 1:N
    ufromx = K_lqr2*(Y3(:,i)-xref(:,i));
    u = uref(:,i) - ufromx;
    u = clipu(u, Fmax, Mmax);
    U3 = [U3,u];
    xnew = Ad2*Y3(:,i)+Bd2*u;
    ynew = xnew;
    Y3 = [Y3,ynew];
end

%%
close all;
hold on;
plot(Y2(1,2:end),Y2(2,2:end));
plot(y(1,:),y(2,:));
title('xref versus simulated');
legend('xsim','xref')
saveas(gcf, 'report/lqr1traj', 'jpg')

figure
hold on;
plot(1:N,U2(1,:))
plot(1:N,U2(2,:))
plot(1:N,U2(3,:))
legend('u1','u2','u3');
title('plot of the input');
saveas(gcf, 'report/lqr1input', 'jpg')

figure
hold on;
plot(1:N,Y2(1,2:end))
plot(1:N,Y2(2,2:end))
plot(1:N,Y2(3,2:end))
plot(1:N,Y2(4,2:end))
plot(1:N,Y2(5,2:end))
plot(1:N,Y2(6,2:end))
legend('x','y','theta','xdot','ydot','thetadot');
saveas(gcf, 'report/lqr1states', 'jpg')

figure;
hold on;
plot(Y3(1,2:end),Y3(2,2:end));
plot(y(1,:),y(2,:));
title('xref versus simulated');
legend('xsim','xref')
saveas(gcf, 'report/lqr2traj', 'jpg')

figure
hold on;
plot(1:N,U3(1,:))
plot(1:N,U3(2,:))
plot(1:N,U3(3,:))
legend('u1','u2','u3');
title('plot of the input');
saveas(gcf, 'report/lqr2input', 'jpg')

figure
hold on;
plot(1:N,Y3(1,2:end))
plot(1:N,Y3(2,2:end))
plot(1:N,Y3(3,2:end))
plot(1:N,Y3(4,2:end))
plot(1:N,Y3(5,2:end))
plot(1:N,Y3(6,2:end))
legend('x','y','theta','xdot','ydot','thetadot');
saveas(gcf, 'report/lqr2states', 'jpg')

%%
plot(0:11,Y2(1,end-1),0:11,xref(1,:))

%% Simulation cost
cost1 = 0;
for i = 1:12
    x = Y2(:,i)-xref(:,i);
    u = U2(:,i);
    cost1 = cost1 + x'*Q1*x +  u'*R1*u;
end

cost2 = 0;
for i = 1:12
    x = Y3(:,i)-xref(:,i);
    u = U3(:,i);
    cost2 = cost2 + x'*Q2*x +  u'*R1*u;
end
cost1
cost2
    
