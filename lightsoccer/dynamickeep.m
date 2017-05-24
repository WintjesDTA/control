%% Initialize system
clear; close all;
m = 70;
ro = 20;
I = 1.4;
Ts = 0.15;
A = [0,0,0,1,0,0;
    0,0,0,0,1,0;
    0,0,0,0,0,1;
    0,0,0,-ro/m,0,0;
    0,0,0,0,-ro/m,0;
    0,0,0,0,0,-.01/I];
B = [0,0,0;
    0,0,0;
    0,0,0;
    1/m,0,0;
    0,1/m,0;
    0,0,1/I];
C = [1,0,0,0,0,0;
    0,1,0,0,0,0;
    0,0,1,0,0,0];
D = zeros(3,3);
sys = ss(A,B,C,D);
sysd = c2d(sys, Ts, 'forward');
step(sysd,0:.15:20)

%% Controllability etc, everything is okay. point 1.2
eig(sysd.A); %Poles
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;
zz=tzero(Ad,Bd,Cd,Dd); %Transmission zeros
%zeta=zz(1); % choose one of the transmission zeros
%M=[zeta*eye(length(Ad))-Ad -Bd; Cd Dd];
%z=null(M);
%x0=z(1:length(Ad),1);
%u0=z(length(Ad)+1:length(z),1);

CC = ctrb(Ad,Bd); %Controllability
[V,dC]=eig(Ad');%unstable uncontrollable -> not stabalizable
Cmodes = Bd'*V;

OO=obsv(Ad,Cd);
rank(OO); % < length(A) ?
[V,dO]=eig(Ad);
Omodes = Cd*V;

sysdeul = sysd;
sysdbilin = c2d(sys, Ts, 'tustin');
sysdzero = c2d(sys, Ts, 'zoh');

inver = inv(eye(6)-A*Ts);
Arect = inver;

Brect = inver*B*Ts;
Crect = C*inver;

Drect = D + C*Brect;
sysdrect = ss(Arect,Brect,Crect,Drect,Ts);
%% simulate systems
tfinal = 6;

subplot(2,2,1)
step(sysdeul,tfinal)
title('Step response of euler');
subplot(2,2,2)
step(sysdbilin,tfinal)
title('Step response of bilinear');
subplot(2,2,3)
step(sysdzero,tfinal)
title('Step response of zero-order');
subplot(2,2,4)
step(sysdrect,tfinal)
title('Step response of rect');
figure;
step(sys,tfinal)
