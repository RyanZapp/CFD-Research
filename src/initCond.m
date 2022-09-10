function [u0,v0] = initCond(x,y,h,J,L,ht)

u0 = zeros(J,J);
v0 = zeros(J,J);
u0(:,1) = -4*y.^3 + 2*y.^2 + 2*y;
v0(:,1) = -(y-0.5).^2 + 0.25;
rev = 1:-h:0;
for j = 2:length(x)
    v0(:,j) = -rev(j)*(y-0.5).^2 + rev(j)*0.25;
    u0(:,j) = rev(j)*(-4*y.^3 + 2*y.^2 + 2*y);
end

