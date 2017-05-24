Fmax = 200;
Mmax = 10;
Nall = size(y,2);
Q = eye(6);
Q(4:6,4:6) = 0;
R = eye(3)*10^-4;
Nref = 12;
[K,S,e] = lqr(sysdeul2,Q,R);

%% Simulate with MPC1
x0 = [0;0;0;0;0;0];
Y4 = [xref(:,1)];
U4 = [];
timings = zeros(Nall,1);
for i = 1:Nall-1
    if i <= Nall-Nref        
        N = Nref;
    else 
        N = Nall-i;
    end
    q = size(R,2);
    p = size(Q,1);
    Qblock = kron(eye(N-1),Q);
    Rblock = kron(eye(N),R);
    a = [Qblock,zeros(p*(N-1),q*N+p)];
    b = [zeros(p,p*(N-1)),S,zeros(p,q*N)];
    c = [zeros(q*N,p*N),Rblock];
    H = 2*[a; b; c];
    Ai = [eye(p*N),zeros(p*N,q*N);
        -eye(p*N),zeros(p*N,q*N);
        zeros(q*N,p*N),eye(q*N);
        zeros(q*N,p*N),-eye(q*N)];
    %xmax = [19.66;16;pi;100;100;100];
    xmax  = [1000000;1000000;1000000;1000000;1000000;1000000];
    umax = [Fmax; Fmax; Mmax];
    bi = [repmat(xmax,2*N,1);repmat(umax,2*N,1)];
    Aepart1 = kron(eye(N),eye(p));
    amatrices = kron(eye(N-1),Ad2);
    Aepart1(p+1:end,1:end-p) = Aepart1(p+1:end,1:end-p)-amatrices;
    Aepart2 = kron(eye(N),-Bd2);
    Ae = [Aepart1,Aepart2];
    Be = [Ad2*Y4(:,i);zeros((N-1)*p,1)];
    temp = [reshape(xref(:,i+1:i+N),N*p,1);reshape(uref(:,i:i+N-1),N*q,1)];
    f = -H*temp;
    tic;
    x = quadprog(H,f,Ai,bi,Ae,Be); 
    timings(i,1) = toc;
    u = x(p*N+1:p*N+3);    
    u = clipu(u, Fmax, Mmax);
    U4 = [U4,u];
    xnew = Ad2*Y4(:,i)+Bd2*u;
    ynew = xnew;
    Y4 = [Y4,ynew];
end

%% Plot stuff
close all;
experiment = 'terminal';
len1 = size(U4,2);
hold on;
plot(Y4(1,2:end),Y4(2,2:end));
plot(y(1,:),y(2,:));
title('xref versus simulated');
legend('xsim','xref')
saveas(gcf, strcat('report/mpctraj',experiment,num2str(Nref)), 'jpg')


figure
hold on;
plot(1:len1,U4(1,:))
plot(1:len1,U4(2,:))
plot(1:len1,U4(3,:))
legend('u1','u2','u3');
title('plot of the input');
saveas(gcf, strcat('report/mpcinput',experiment,num2str(Nref)), 'jpg')

figure
hold on;
plot(1:len1,Y4(1,2:end))
plot(1:len1,Y4(2,2:end))
plot(1:len1,Y4(3,2:end))
plot(1:len1,Y4(4,2:end))
plot(1:len1,Y4(5,2:end))
plot(1:len1,Y4(6,2:end))
legend('x','y','theta','xdot','ydot','thetadot');
saveas(gcf, strcat('report/mpcstates',experiment,num2str(Nref)), 'jpg')
%%
simcost = 0;
for i = 1:Nall-1
    x = Y4(:,i)-xref(:,i);
    u = U4(:,i);
    simcost = simcost + x'*Q2*x +  u'*R1*u;
end
simcost
%%
figure;
plot(1:size(timings,1),timings)
xlabel('iteration number')
ylabel('computation time in sec')
title('computational complexity plot')
axis([0 160 0 .03])
saveas(gcf, strcat('report/comptime',experiment,num2str(Nref)), 'jpg')
%%
xmax = [19.66;19.66;pi;100;100;100];
umax = [Fmax; Fmax; Mmax];
bi = [repmat(xmax,2*N,1);repmat(umax,2*N,1)];
Aepart1 = kron(eye(N),eye(p));
amatrices = kron(eye(N-1),Ad2);
Aepart1(p+1:end,1:end-p) = -amatrices;
Aepart2 = kron(eye(N),-Bd2);
Ae = [Aepart1,Aepart2];