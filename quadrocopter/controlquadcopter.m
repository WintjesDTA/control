%%
load('references_17.mat')



m = .5;
L = .25;
k = 3*10^-6;
b = 10^-7;
g = 9.81;
kd = .25;
cm = 10^4;
Ixx = .005;
Iyy = .005;
Izz = .01;

hv = g*m/(4*k*cm); % hover_voltage; 

A = [0,0,0,1,zeros(1,8);
    0,0,0,0,1,zeros(1,7);
    0,0,0,0,0,1,zeros(1,6);
    0,0,0,-kd/m,zeros(1,8);
    0,0,0,0,-kd/m,zeros(1,7);
    0,0,0,0,0,-kd/m,zeros(1,6);
    zeros(1,9),1,0,0;
    zeros(1,10),1,0;
    zeros(1,11),1;
    zeros(1,12);
    zeros(1,12);
    zeros(1,12)];  
A(5,7) = -g;
A(4,8) = g;
B = zeros(12,4);
B(6,:) = (k*cm/m)*ones(4,1); %vx input
B(10,1) = L*k*cm/Ixx; B(10,3) = -L*k*cm/Ixx;
B(11,2) = L*k*cm/Iyy; B(11,4) = -L*k*cm/Iyy;
B(12,[1,3]) = b*cm/Izz; B(12,[2,4]) = -b*cm/Izz;
%B(4,:) = ones(4,1); %vy input
C = zeros(6,12);
C(1:3,1:3) = eye(3);
C(4:6,7:9) = eye(3);
D = zeros(6,4);
sysl = ss(A,B,C,D);
t = linspace(0,5,1000);
u = zeros(4,length(t));
u(:,200:1000) = u(:,200:1000) + 5;
lsim(sysl,u,t);

%% Discrete system, stability, poles, controllability.
Ts = 0.05;
sysd = c2d(sysl,Ts); 

eig(sysd.A) %Poles
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;
zz=tzero(Ad,Bd,Cd,Dd); %Transmission zeros
zeta=zz(1); % choose one of the transmission zeros
M=[zeta*eye(length(Ad))-Ad -Bd; Cd Dd];
z=null(M);
x0=z(1:length(Ad),1);
u0=z(length(Ad)+1:length(z),1);

CC = ctrb(Ad,Bd); %Controllability
[V,dC]=eig(Ad');%unstable uncontrollable -> not stabalizable
Cmodes = Bd'*V;

OO=obsv(Ad,Cd);
rank(OO); % < length(A) ?
[V,dO]=eig(Ad);
Omodes = Cd*V;

%% LQR control
Ts = 0.05;
%Q = zeros(12,12);
%Q(1:3,1:3) = [1,0,0; 0,1,0;0,0,1];
%Q = 1000*eye(12);
%Q(1:3,1:3) = 10000*eye(3);
%R = 1*eye(4);
Q = 100*eye(12);
Q(1:2,1:2) = 500*eye(2); %500
Q(3,3) = 50000;
%Q(1:3,1:3) = 1000*eye(3); 
%Q(3,3) = 10000;
R = 1*eye(4);
a = 1;
b = 1;
Cnew = C(1:3,:); Dnew = D(1:3,:);
sysd = c2d(ss(A,B,Cnew,Dnew),Ts, 'zoh'); %Was tustin first
K_lqr = lqr(sysd,Q,R,0); 
Z = [sysd.A-eye(12),sysd.B;sysd.C,sysd.D];
N = pinv(Z)*[zeros(12,3);eye(3)];
Nx = N(1:12,:);
Nu = N(13:16,:);
Nbar = Nu + K_lqr*Nx;
K_lqr

%% Integral control execute LQR first
Q = 1500*eye(15);
%Q(1:3,1:3) = 500;
Q(3,3) = 3000;
%Q(1:3,1:3) = 10000*eye(3);
Q(13:15,13:15) = 40000*eye(3); %30000 no payload
Q(13:14,13:13) = 50000;
R = 1*eye(4);
[K_lqi,S,e] = lqi(sysd,Q,R,0);
K_lqi

%%
sysd2 = c2d(ss(A,B,C,D),Ts, 'zoh');
Q_kal = zeros(12,12);
Q_kal(4:6,4:6) = eye(3);
Q_kal(10:12,10:12) = eye(3);
Q_kal = Q_kal*0.01;
R_kal = zeros(6,6);
R_kal(1:3,1:3) = 2.5*10^-5*eye(3);
R_kal(4:6,4:6) = 7.57*10^-5*eye(3);