K = 1;
tf = 3;
Nt = 100;
t = linspace(0,tf,Nt);
dt = (tf-0)/(Nt-1);
a = 0;
b = pi;
Nx = 1000;
x = linspace(a,b,Nx);
dx = (b-a)/(Nx-1);
% set initial condition
u = sin(x);
% set boundary conditions
u(1) = 0;
u(end) = 0;
% set test functions (we will use 2 test functions
phi1 = (x-a)/(b-a);
phi1 = sin(x);
phi2 = (x-b)/(a-b);
% compute the derivatives of the test functions
dphi1dx = diff(phi1);
dphi1dx(end+1) = dphi1dx(end);
dphi2dx = diff(phi2);
dphi2dx(end+1) = dphi2dx(end);
G11_vect = (phi1.*phi1 + K*dt*dphi1dx.*dphi1dx)*dx;
G12_vect = (phi1.*phi2 + K*dt*dphi1dx.*dphi2dx)*dx;
G21_vect = (phi2.*phi1 + K*dt*dphi2dx.*dphi1dx)*dx;
G22_vect = (phi2.*phi2 + K*dt*dphi2dx.*dphi2dx)*dx;
G11 = sum(G11_vect);
G12 = sum(G12_vect);
G21 = sum(G21_vect);
G22 = sum(G22_vect);
G = [G11 G12; G21 G22];
G_inv = inv(G);
figure
count = 1;
for i = 1:Nt
    F1_vect = u.*phi1*dx;
    F2_vect = u.*phi2*dx;
    F1 = sum(F1_vect);
    F2 = sum(F2_vect);
    F = [F1; F2];
    % compute the vector and call it F
    C = G_inv*F;
    u = C(1)*phi1 + C(2)*phi2; % THis is the new u value at t+1
    plot(x,u)
    axis([0 pi -inf inf])
    drawnow
    count = count+1
end
% Its definitely not legit because it
%figure
%plot(x,phi1,'b')
%hold on
%plot(x,phi2,'r')
%axis([0 pi -inf inf])