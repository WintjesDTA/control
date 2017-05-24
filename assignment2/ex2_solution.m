clear
close all

load ex2_data

[n_states,n_inputs] = size(B2);
% n_outputs = size(C,1);

%% question 1

G = randn(n_inputs,n_states);

Lambda1 = blkdiag([-3 +1; -1 -3],-1);
X1 = lyap(A,-Lambda1,-B2*G);
K1 = G/X1
% check
eig(A-B2*K1)

Lambda2 = blkdiag([-10 +2.5; -2.5 -10],-12.1);
X2 = lyap(A,-Lambda2,-B2*G);
K2 = G/X2
% check
eig(A-B2*K2)

sys1 = ss(A-B2*K1,B1(:,1),[C2;C1;-D1*K1;-K1],0);
figure
step(sys1,10);
sys2 = ss(A-B2*K2,B1(:,1),[C2;C1;-D1*K2;-K2],0);
figure
step(sys2,10);

pause

%% question 2

P = are(A,B2/(D1'*D1)*B2',C1'*C1);
K = (D1'*D1)\B2'*P

K_lqr = lqr(A,B2,C1'*C1,D1'*D1)

Q1 = 10^(-3)*eye(n_states);
Q2 = 10^(3)*eye(n_states);
R  = eye(n_inputs);
K1 = lqr(A,B2,Q1,R);
K2 = lqr(A,B2,Q2,R);
sys1 = ss(A-B2*K1,B1(:,1),[C2;C1;-D1*K1;-K1],0);
figure
impulse(sys1,10);
sys2 = ss(A-B2*K2,B1(:,1),[C2;C1;-D1*K2;-K2],0);
figure
impulse(sys2,10);

pause

%% question 3

L1 = place(A',C2',eig(A-B2*K_lqr)*1.1)'
L2 = place(A',C2',eig(A-B2*K_lqr)*5)'

T = 999;
w = 0.1*randn(size(B1,2),T);
v = 0.1*randn(size(L1,2),T);

sys1 = ss(A-L1*C2,[B1 -L1],eye(size(A)),zeros(size(A,1),size(B1,2)+size(L1,2)));
t = (0:1:(T-1))/100;
error1 = lsim(sys1,[w;v],t,[1;1;1]);
sys2 = ss(A-L2*C2,[B1 -L2],eye(size(A)),zeros(size(A,1),size(B1,2)+size(L2,2)));
t = (0:1:(T-1))/100;
error2 = lsim(sys2,[w;v],t,[1;1;1]);

sys_x = ss(A,zeros(n_states,1),eye(n_states,n_states),zeros(n_states,1));
x = lsim(sys_x,zeros(size(t)),t,[1;1;1]);

figure
for i=1:3
    subplot(3,2,2*i-1)
    plot(t,[x(:,i) x(:,i)-error1(:,i)])
    legend('x','x est')
    subplot(3,2,2*i)
    plot(t,error1(:,i));
end
figure
for i=1:3
    subplot(3,2,2*i-1)
    plot(t,[x(:,i) x(:,i)-error2(:,i)])
    legend('x','x est')
    subplot(3,2,2*i)
    plot(t,error2(:,i));
end

pause

%% question 4

L = lqe(A,B1,C2,eye(size(B1,2)),eye(size(C2,1)))

L1 = lqe(A,B1,C2,0.1*eye(size(B1,2)),0.1*eye(size(C2,1)))
L2 = lqe(A,B1,C2,0.01*eye(size(B1,2)),0.1*eye(size(C2,1)))
sys1 = ss(A-L1*C2,[B1 -L1],eye(size(A)),zeros(size(A,1),size(B1,2)+size(L1,2)));
t = (0:1:(T-1))/100;
figure
lsim(sys1,[w;v],t,[1;1;1])
sys2 = ss(A-L2*C2,[B1 -L2],eye(size(A)),zeros(size(A,1),size(B1,2)+size(L2,2)));
t = (0:1:(T-1))/100;
figure
lsim(sys2,[w;v],t,[1;1;1])

T = 999;
w = 1*randn(size(B1,2),T);
v = 0.1*randn(size(L1,2),T);

L1 = lqe(A,B1,C2,1*eye(size(B1,2)),0.1*eye(size(C2,1)))
L2 = lqe(A,B1,C2,0.1*eye(size(B1,2)),1*eye(size(C2,1)))
sys1 = ss(A-L1*C2,[B1 -L1],eye(size(A)),zeros(size(A,1),size(B1,2)+size(L1,2)));
t = (0:1:(T-1))/100;
error1 = lsim(sys1,[w;v],t,[1;1;1]);
sys2 = ss(A-L2*C2,[B1 -L2],eye(size(A)),zeros(size(A,1),size(B1,2)+size(L2,2)));
t = (0:1:(T-1))/100;
error2 = lsim(sys2,[w;v],t,[1;1;1]);

sys_x = ss(A,zeros(n_states,1),eye(n_states,n_states),zeros(n_states,1));
x = lsim(sys_x,zeros(size(t)),t,[1;1;1]);

figure
for i=1:3
    subplot(3,2,2*i-1)
    plot(t,[x(:,i) x(:,i)-error1(:,i)])
    legend('x','x est')
    subplot(3,2,2*i)
    plot(t,error1(:,i));
end
figure
for i=1:3
    subplot(3,2,2*i-1)
    plot(t,[x(:,i) x(:,i)-error2(:,i)])
    legend('x','x est')
    subplot(3,2,2*i)
    plot(t,error2(:,i));
end