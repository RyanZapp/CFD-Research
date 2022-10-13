clear
close all
clc
%pp
L = 2;
ht = 1;
rho = 1 ;

t_final = 2 ;
N = 10000 ;
J = 250 ;

nu = 0.008926 ;
mu = rho*nu ;
% Eventually swap h with dx and dy
h = ht/(J-1) ;
deltaT = 0.125/nu*h^2 ;

x = linspace(0,L,J) ;
y = linspace(0,ht,J) ;
k = -0.5*x + 1;
t = linspace(0,t_final,N);

u = zeros(J,J,N);
v = zeros(J,J,N);

[u0,v0] = initCond(x,y,h,J,L,ht);

for i = 2:N
    u(:,1,i) = u0(:,1);
    v(:,1,i) = v0(:,1);
end