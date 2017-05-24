function [u] = clipu(u,Fmax,Mmax)
%CLIPU Summary of this function goes here
%   Detailed explanation goes here
    if abs(u(1)) > Fmax
        u(1) = Fmax*sign(u(1));           
    end
    if abs(u(2)) > Fmax
        u(2) = Fmax*sign(u(2));           
    end
    if abs(u(3)) > Mmax
        u(3) = 10*sign(u(3));
    end
end

