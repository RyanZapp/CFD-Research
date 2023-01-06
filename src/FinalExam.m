%%Problem 1 Volume Integral
N = 200;
[x,y,z] = meshgrid(linspace(0,1,N));
dX = (1)/(N-1);
dY = (1)/(N-1);
dZ = (1)/(N-1);
dV = dX*dY*dZ;
f = x;
g = y;
h = z;
div = divergence(x,y,z,f,g,h);
Flux = sum(div(:))*dV;
disp('Flux from volume integral:')
disp(Flux)
%% Problem 1 Surface Integral
clear
close all
clc
N = 5;
[x,y,z] = meshgrid(linspace(0,1,5));
f = x;
g = y;
m = z;
[u,v] = meshgrid(linspace(0,1,N));
J = size(u);
oo = zeros(J(1),J(2));
k = ones(J(1),J(2));
x1 = cat(3,k,u,v);
x2 = cat(3,u,k,v);
x3 = cat(3,oo,u,v);
x4 = cat(3,u,oo,v);
x5 = cat(3,u,v,oo);
x6 = cat(3,u,v,k);
figure
h = surf(x1(:,:,1),x1(:,:,2),x1(:,:,3));
set(h,'FaceColor','b')
hold on
h = surf(x2(:,:,1),x2(:,:,2),x2(:,:,3));
set(h,'FaceColor','b')
hold on
h = surf(x3(:,:,1),x3(:,:,2),x3(:,:,3));
set(h,'FaceColor','b')
hold on
h = surf(x4(:,:,1),x4(:,:,2),x4(:,:,3));
set(h,'FaceColor','b')
hold on
h = surf(x5(:,:,1),x5(:,:,2),x5(:,:,3));
set(h,'FaceColor','b')
hold on
h = surf(x6(:,:,1),x6(:,:,2),x6(:,:,3));
set(h,'FaceColor','b')
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Plot of Field and Domain')

a1 = diff(x1,1,1);
a1(:,end,:) = [];
b1 = diff(x1,1,2);
b1(end,:,:) = [];
n1 = cross(a1,b1,3);
dA1 = sqrt(sum(n1.^2,3));

a2 = diff(x2,1,1);
a2(:,end,:) = [];
b2 = diff(x2,1,2);
b2(end,:,:) = [];
n2 = cross(a2,b2,3);
dA2 = sqrt(sum(n2.^2,3));

a3 = diff(x3,1,1);
a3(:,end,:) = [];
b3 = diff(x3,1,2);
b3(end,:,:) = [];
n3 = cross(a3,b3,3);
dA3 = sqrt(sum(n3.^2,3));

a4 = diff(x4,1,1);
a4(:,end,:) = [];
b4 = diff(x4,1,2);
b4(end,:,:) = [];
n4 = cross(a4,b4,3);
dA4 = sqrt(sum(n4.^2,3));

a5 = diff(x5,1,1);
a5(:,end,:) = [];
b5 = diff(x5,1,2);
b5(end,:,:) = [];
n5 = cross(a5,b5,3);
dA5 = sqrt(sum(n5.^2,3));

a6 = diff(x6,1,1);
a6(:,end,:) = [];
b6 = diff(x6,1,2);
b6(end,:,:) = [];
n6 = cross(a6,b6,3);
dA6 = sqrt(sum(n6.^2,3));

for j = 1:3
    n1(:,:,j) = -n1(:,:,j)./dA1;
end

for j = 1:3
    n2(:,:,j) = n2(:,:,j)./dA2;
end

for j = 1:3
    n3(:,:,j) = n3(:,:,j)./dA3;
end

for j = 1:3
    n4(:,:,j) = -1*n4(:,:,j)./dA4;
end

for j = 1:3
    n5(:,:,j) = n5(:,:,j)./dA5;
end

for j = 1:3
    n6(:,:,j) = -n6(:,:,j)./dA6;
end

x1(end,:,:) = [];
x1(:,end,:) = [];

x2(end,:,:) = [];
x2(:,end,:) = [];

x3(end,:,:) = [];
x3(:,end,:) = [];

x4(end,:,:) = [];
x4(:,end,:) = [];

x5(end,:,:) = [];
x5(:,end,:) = [];

x6(end,:,:) = [];
x6(:,end,:) = [];

x(end,:,:) = [];
x(:,end,:) = [];
x(:,:,(end-1):end) = [];

y(end,:,:) = [];
y(:,end,:) = [];
y(:,:,(end-1):end) = [];

z(end,:,:) = [];
z(:,end,:) = [];
z(:,:,(end-1):end) = [];


f(end,:,:) = [];
f(:,end,:) = [];
f(:,:,(end-1):end) = [];

g(end,:,:) = [];
g(:,end,:) = [];
g(:,:,(end-1):end) = [];

m(end,:,:) = [];
m(:,end,:) = [];
m(:,:,(end-1):end) = [];

quiver3(x1(:,:,1),x1(:,:,2),x1(:,:,3),n1(:,:,1),n1(:,:,2),n1(:,:,3),0,'k')
hold on
quiver3(x2(:,:,1),x2(:,:,2),x2(:,:,3),n2(:,:,1),n2(:,:,2),n2(:,:,3),0,'k')
hold on
quiver3(x3(:,:,1),x3(:,:,2),x3(:,:,3),n3(:,:,1),n3(:,:,2),n3(:,:,3),0,'k')
hold on
quiver3(x4(:,:,1),x4(:,:,2),x4(:,:,3),n4(:,:,1),n4(:,:,2),n4(:,:,3),0,'k')
hold on
quiver3(x5(:,:,1),x5(:,:,2),x5(:,:,3),n5(:,:,1),n5(:,:,2),n5(:,:,3),0,'k')
hold on
quiver3(x6(:,:,1),x6(:,:,2),x6(:,:,3),n6(:,:,1),n6(:,:,2),n6(:,:,3),0,'k')
hold on
quiver3(x,y,z,f,g,m,0,'r')
hold on

f1 = x1;
f2 = x2;
f3 = x3;
f4 = x4;
f5 = x5;
f6 = x6;

flux1 = sum(f1.*n1,3).*dA1;
flux1 = sum(flux1(:));

flux2 = sum(f2.*n2,3).*dA2;
flux2 = sum(flux2(:));

flux3 = sum(f3.*n3,3).*dA3;
flux3 = sum(flux3(:));

flux4 = sum(f4.*n4,3).*dA4;
flux4 = sum(flux4(:));

flux5 = sum(f5.*n5,3).*dA5;
flux5 = sum(flux5(:));

flux6 = sum(f6.*n6,3).*dA6;
flux6 = sum(flux6(:));

Flux = flux1 + flux2 + flux3 + flux4 + flux5 + flux6;
disp('Flux computed via surface integral:')
disp(Flux)
%% Problem 2

clear
close all
clc
% Parameterize Triangle
eps = 0.1;
N = 2500; % Increase later
t = linspace(0,1,N)';
A = 0.5*eps^2;

x1 = eps*t;
y1 = 0*t; % multiply by t to make y1 the same size as x1
x2 = eps*(1-t);
y2 = eps*t;
x3 = 0*t;
y3 = eps*(1-t);

f1 = [-y1, x1];
f2 = [-y2, x2];
f3 = [-y3, x3];

% Compute tangent vector
% Then append final tangent vector to end of diff
dx1 = diff(x1); % computes forward difference
dy1 = diff(y1);
% We can set Nth tangent vector = N-1th tangent vector
dx1(end+1) = dx1(end);
dy1(end+1) = dy1(end);
dx2 = diff(x2);
dy2 = diff(y2);
dx2(end+1) = dx2(end);
dy2(end+1) = dy2(end);
dx3 = diff(x3);
dy3 = diff(y3);
dx3(end+1) = dx3(end);
dy3(end+1) = dy3(end);

% Store components of tangent vector for each side of triangle
ds1 = [dx1,dy1];
ds2 = [dx2,dy2];
ds3 = [dx3,dy3];
% Compute unit tangent vector for each side
t1 = ds1./vecnorm(ds1,2,2);
t2 = ds2./vecnorm(ds2,2,2);
t3 = ds3./vecnorm(ds3,2,2);
% Store and sum f dot t *ds for each segment along the boundary
f1_dot_t1_ds = f1.*t1.*vecnorm(ds1,2,2);
f2_dot_t2_ds = f2.*t2.*vecnorm(ds2,2,2);
f3_dot_t3_ds = f3.*t3.*vecnorm(ds3,2,2);
work1 = sum(f1_dot_t1_ds(:));
work2 = sum(f2_dot_t2_ds(:));
work3 = sum(f3_dot_t3_ds(:));
work = 1/A*(work1 + work2 + work3);
disp('Work divided by area from boundaryintegral:')
disp(work)
figure
plot(x1,y1)
hold on
plot(x2,y2)
hold on
plot(x3,y3)