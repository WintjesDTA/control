[X,Y] = meshgrid(-30:.5:30,-10:.5:50);
num = size(X,1);
Z = zeros(num,num);
for i = 1:num
    for j = 1:num
        Z(i,j) = costfunction(X(i,j),Y(i,j),0);
    end
end
Z = log(Z);
surf(X,Y,Z)
xlabel('x-coord');
ylabel('y-coord');  
colorbar
title('Plot of the costfunction');