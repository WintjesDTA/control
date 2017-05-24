syms ro m I Ts
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
Adeul = eye(6) + A*Ts;
Bdeul = B*Ts;
Cdeul = C;
Ddeul = D;
%% Fill in numerical values
m = 70;
ro = 20;
I = 1.4;
Ts = 0.15;
Ad = double(subs(Adeul));
Bd = double(subs(Bdeul));
Cd = Cdeul;
Dd = Ddeul;
Ts = double(subs(Ts));
%sysd = c2d(sys, Ts, 'forward');