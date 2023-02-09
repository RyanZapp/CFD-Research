clear
close all
clc
N = 100;
%B = [-2,0,1,4];
y = linspace(-1,1,N);
%u = 0.5*(-B*y.^2 + y + B + 1);
figure
for B = [-2,0,1,4]
    u = 0.5*(-B*y.^2 + y + B + 1);
    plot(u,y)
    hold on
end
grid on
title('Velocity profile for different values of B')
legend('B=-2','B=0','B=1','B=4','location','SE')
xlabel('u*')
ylabel('y*')
% I should probably mark the y = -1/2 and y = 1/2 points to make plotting
% easier (Or I can ask Lutz if it is ok to do it in MATLAB)
axis([-1 3 -1 1])