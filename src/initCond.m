function [u0,v0] = initCond(x,y,dxi,dyi,nx,ny,Lx,Ly)

u0 = zeros(nx,ny);
v0 = zeros(nx,ny);
u0(:,1) = -4*y.^3 + 2*y.^2 + 2*y; % What am I going to do in order to properly set the initial and boundary conditions
v0(:,1) = -(y-0.5).^2 + 0.25;
rev = 1:-dyi:0;
for j = 2:length(x)
    v0(:,j) = -rev(j)*(y-0.5).^2 + rev(j)*0.25;
    u0(:,j) = rev(j)*(-4*y.^3 + 2*y.^2 + 2*y);
end

