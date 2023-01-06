clear
close all
clc
format
% Before reading in the file, delete the name of the airfoil from the file
% header
T = readtable('NACA0024.txt');
A = table2array(T);
x = A(:,1);
y = A(:,2);
xo = x;
yo = y;
xq = 1; % finer querey range
yq = 1; % finer querey range

figure
plot(x,y)
axis equal