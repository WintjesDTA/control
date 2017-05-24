minx = -16 - 7.32/2;
maxx = 16 + 7.32/2;
LB = [minx,0,-pi/2];
UB = [maxx,16,pi/2];
A = zeros(3,3);
b = zeros(3,1);
x0  = [1,1,1];
X = fmincon(@costfunction,x0,A,b,A,b,LB,UB);

