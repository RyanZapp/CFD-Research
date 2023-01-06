clear
close all
clc
Nx = 5;
Ny = 5;
Lx = 1;
Ly = 1;
dx = Lx/Nx;
dy = Ly/Ny; % This is basically (b-a)/N
xce = ((1:Nx)-0.5)*dx; % this is basically doing x_j = a + dx*j
% but in this case, a = 0 and we are shifting x_j back by 0.5 as well
yce = ((1:Ny)-0.5)*dy;
u = zeros(Nx,Ny);
u(:,1) = 0; % after the transform, this should be the bottom side
u(:,end) = 1; % after the transform, this should be the top
u(1,:) = 3; % after the transform, this should be the left
u(end,:) = 5; % after the transform, this should be the right
[Xce,Yce] = meshgrid(xce,yce);
X = Xce';
Y = Yce';

% steal the code below to also try and figure out what is good, by
% transposing things in the contourf plot
% x = linspace(-2*pi,2*pi);
% y = linspace(0,4*pi);
% [X,Y] = meshgrid(x,y);
% Z = sin(X) + cos(Y);
% contourf(X,Y,Z,10)