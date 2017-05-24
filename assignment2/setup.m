load('ex2_data.mat');
lambda = eig(A); %Get the poles of the system.

%% find the matrix Gamma
alpha = real(lambda(2));
beta = imag(lambda(2));
lambda1 = real(lambda(1));
gamma = [alpha, beta, 0;
        -beta, alpha, 0;
        0, 0, lambda1];
K1 = place(A,B2,[-3+1i,-3-1i, -1]);
K2 = place(A,B2,[-10+2.5*1i,-10-2.5*1i, -12.1]);
eig(A-B2*K1) %works
eig(A-B2*K2) %works

%% Get impulse responses
hold on;
sys  = ss(A,B2,C2,0);
sys1 = ss(A-B2*K1,B2,C2,0);
sys2 = ss(A-B2*K2,B2,C2,0);
impulse(sys1,sys2);
%impulse(sys2);
legend('system1','system2')

%how to add noise to the system, check solution.

%% 2. Control law design
n = .2;
B = B2*B2' + n*eye(3);
C = C1'*C1 + n*eye(3);
P = are(A,B,C);
%K = eye(2)*B'*P;

K_lqr = lqr(A, B2, C, n*eye(2));
syslqr = ss(A-B2*K_lqr,B2,C2,0);
impulse(syslqr)
eig(A-B2*K_lqr)

%% 3. Observer Design
% a)
L1 = place(A',C2',eig(A-B2*K_lqr)*1.1)';
L2 = place(A',C2',eig(A-B2*K_lqr)*5)';

% b) create white noise
w = randn(300);
v = randn(300);
t = 0:.01:3;
xtil0 = [0,0,0];

Aobs = (A-L*C2);
Bobs = [B1,-L];
sysobs = ss(Aobs, Bobs, eye(2), 0);
y = lsim(sysobs,[w;v],t,xtil0);





