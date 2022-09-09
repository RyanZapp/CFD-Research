clear
close all
clc
I = 100;
J = 100;
N = 100; % Place the function so that it can read in I and J from the main matlab file and use those to create the size of the initial condition matrix
% Really you should be trying to interface this so that it reads in values
% from the navier stokes code
k = 1;
y = linspace(0,1,J);
f = -k*(y-0.5).^2 +0.25*k; % This is our boundary condition for the inlet I believe
v0 = zeros(I,J);
for j = 1:J % fix this so that the value of the y function decreases over time
    v0(:,j) = -abs(j)*(y-0.5).^2 + 0.25*abs(j);
end
%figure
%plot(f,y)