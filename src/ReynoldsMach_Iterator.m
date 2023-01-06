clear
close all
clc
% All values are given in m-kg-s
% Temperature dependent constants are currently being evaluated at 60
% degrees F
% format long g
format
tsr = 2;
mu = 17.91*10^-6; % N*s/(m^2)
rho = 1.222; % kg/(m^3)
Chord_Length = 0.0762; % (m)
a = 344; % (m/s) Speed of sound in air
N = 5; % Number of velocity sample points
v_min = 33.528; % (m/s)
v_max = tsr*58.1152; % (m/s)
delta_v = (v_max - v_min)/(N-1);
delta_Re = rho*delta_v*Chord_Length/mu; % Change in the Reynolds number
v = linspace(v_min,v_max,N);
Re = rho*v*Chord_Length/mu; % Reynolds Number (Dimensionless)
M = v/a; % Mach Number (Dimensionless)

fid1 = 'Reynolds Number Lower Bound: %.5f \n';
fprintf(fid1,Re(1))
fid2 = 'Reynolds Number Upper Bound: %.5f \n';
fprintf(fid2,Re(end))
fid3 = 'Change in Reynolds Number: %.5f \n';
fprintf(fid3,delta_Re)
fprintf('\n');
fid4 = 'Mach Number Upper Bound: %.5f \n';
fprintf(fid4,M(end))

