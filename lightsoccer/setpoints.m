%% Set up input signal, run dynamic keep first
load('trajectory_gk.mat')
N = size(y,2);
sysdopt = sysdeul;
Ad = sysdopt.A;
Bd = sysdopt.B;
Cd = sysdopt.C;
Dd = sysdopt.D;
H = zeros(N,3,3);

ON = [];
for k = 0:N-1
    ON = [ON;Cd*Ad^k];
end

len1 = size(D,1);
len2 = size(D,2);
HN = zeros(N*len1,N*len2);
for k = 0:N-1
    if k == 0
        temp = Dd;
    else
        temp = Cd*Ad^(k-1)*Bd;
    end
    addition = kron(eye(N-k),temp);
    index1 = k*len1+1;
    index2 = (N-k)*len2;    
    HN(index1:end,1:index2) = HN(k*len1+1:end,1:(N-k)*len2) + addition;
end

% D = 0, we use quadprog

x0 = [0;0;0;0;0;0];
w = [];
for i = 1:N
    w = [w;y(1,i);y(2,i);y(3,i)];
end
mu= 0.000001;  %0.00001 much better
H = 2*HN'*HN + mu*eye(size(HN'*HN));
f = -2*w'*HN + 2*x0'*ON'*HN;
[se,val] = quadprog(H,f,zeros(3*N,3*N),zeros(3*N,1));

uopt = [];
for i = 0:N-1
    uopt = [uopt,[se(i*len1+1);se(i*len1+2);se(i*len1+3)]];
end

Y = [x0];
for i = 1:N
    xnew = Ad*Y(:,i)+Bd*uopt(:,i);
    ynew = xnew;
    Y = [Y,ynew];
end


%% Plot signals
close all;
figure;
hold on;
plot(Y(1,2:end),Y(2,2:end));
plot(y(1,:),y(2,:));
title('xref versus simulated');
legend('xsim','xref')
saveas(gcf, 'report/xrefstates','jpg');

figure
hold on;
plot(1:N,uopt(1,:))
plot(1:N,uopt(2,:))
plot(1:N,uopt(3,:))
legend('u1','u2','u3');
title('plot of the input');
saveas(gcf, 'report/inputwithreg','jpg');

%%
plot(0:11,Y(1,2:end),0:11,Y(2,2:end),0:11,Y(3,2:end))
legend('xref','yref','thetaref');





