clear
close all
clc
% Shrinkable parameters
N = 10000; % makes parameterization more accurate
dr = 0.000001; % shrinking this is analogous to taking the limit as dr -> 0
dtheta = 0.000001; % shrinking this is analogous to taking the limit as dtheta -> 0

b = 0.5;
a = -0.5;
t = linspace(a,b,N);
dt = (b-a)/(N-1); % we can compute dt rather than using the diff command on t
dti = 1/dt;

theta0 = pi/4;
r0 = 5;

dA = r0*dr*dtheta;



r_hat = [cos(theta0); sin(theta0)];
theta_hat = [-sin(theta0); cos(theta0)];



s_T_r = (r0+0.5*dr)*cos(dtheta*t).*r_hat;
s_T_theta = (r0+0.5*dr)*sin(dtheta*t).*theta_hat;
s_T = s_T_r + s_T_theta;

s_B_r = (r0-0.5*dr)*cos(dtheta*t).*r_hat;
s_B_theta = -(r0-0.5*dr)*sin(dtheta*t).*theta_hat;
s_B = s_B_r + s_B_theta;

s_R_r = (r0+dr*t)*cos(dtheta/2).*r_hat;
s_R_theta = -(r0+dr*t)*sin(dtheta/2).*theta_hat;
s_R = s_R_r + s_R_theta;


s_L_r = (r0-dr*t)*cos(dtheta/2).*r_hat;
s_L_theta = (r0-dr*t)*sin(dtheta/2).*theta_hat;
s_L = s_L_r + s_L_theta;
% first row is the x component and 2nd row is the y component

% my parameterizations look quite good

% Compute analytical derivatives:
ds_T_rdt = -dtheta*(r0+0.5*dr)*sin(dtheta*t).*r_hat;
ds_T_thetadt = dtheta*(r0+0.5*dr)*cos(dtheta*t).*theta_hat;
ds_T_dt = ds_T_rdt + ds_T_thetadt;

ds_B_rdt = -dtheta*(r0-0.5*dr)*sin(dtheta*t).*r_hat;
ds_B_thetadt = -dtheta*(r0-0.5*dr)*cos(dtheta*t).*theta_hat;
ds_B_dt = ds_B_rdt + ds_B_thetadt;

ds_R_rdt = dr*cos(0.5*dtheta).*r_hat;
ds_R_thetadt = -dr*sin(0.5*dtheta).*theta_hat;
ds_R_dt = ds_R_rdt + ds_R_thetadt;
for i = 1:N
    ds_R_dt(:,i) = ds_R_dt(:,1);
end

ds_L_rdt = -dr*cos(0.5*dtheta).*r_hat;
ds_L_thetadt = -dr*sin(0.5*dtheta).*theta_hat;
ds_L_dt = ds_L_rdt + ds_L_thetadt;
for i = 1:N
    ds_L_dt(:,i) = ds_L_dt(:,1);
end

% Compute numerical derivatives
dsdt_T_r = [diff(s_T_r(1,:)); diff(s_T_r(2,:))];
dsdt_T_r(:,end+1) = [dsdt_T_r(1,end); dsdt_T_r(2,end)];
dsdt_T_r = dsdt_T_r*dti;

dsdt_T_theta = [diff(s_T_theta(1,:)); diff(s_T_theta(2,:))];
dsdt_T_theta(:,end+1) = [dsdt_T_theta(1,end); dsdt_T_theta(2,end)];
dsdt_T_theta = dsdt_T_theta*dti;

dsdt_T = dsdt_T_r + dsdt_T_theta; % numerical top derivative


dsdt_B_r = [diff(s_B_r(1,:)); diff(s_B_r(2,:))];
dsdt_B_r(:,end+1) = [dsdt_B_r(1,end); dsdt_B_r(2,end)];
dsdt_B_r = dsdt_B_r*dti;

dsdt_B_theta = [diff(s_B_theta(1,:)); diff(s_B_theta(2,:))];
dsdt_B_theta(:,end+1) = [dsdt_B_theta(1,end); dsdt_B_theta(2,end)];
dsdt_B_theta = dsdt_B_theta*dti;

dsdt_B = dsdt_B_r + dsdt_B_theta; % numerical bottom derivative

dsdt_R_r = [diff(s_R_r(1,:)); diff(s_R_r(2,:))];
dsdt_R_r(:,end+1) = [dsdt_R_r(1,end); dsdt_R_r(2,end)];
dsdt_R_r = dsdt_R_r*dti;

dsdt_R_theta = [diff(s_R_theta(1,:)); diff(s_R_theta(2,:))];
dsdt_R_theta(:,end+1) = [dsdt_R_theta(1,end); dsdt_R_theta(2,end)];
dsdt_R_theta = dsdt_R_theta*dti;

dsdt_R = dsdt_R_r + dsdt_R_theta; % numerical right derivative
dsdt_R_Plot = [mean(dsdt_R(1,:)); mean(dsdt_R(2,:))];

dsdt_L_r = [diff(s_L_r(1,:)); diff(s_L_r(2,:))];
dsdt_L_r(:,end+1) = [dsdt_L_r(1,end); dsdt_L_r(2,end)];
dsdt_L_r = dsdt_L_r*dti;

dsdt_L_theta = [diff(s_L_theta(1,:)); diff(s_L_theta(2,:))];
dsdt_L_theta(:,end+1) = [dsdt_L_theta(1,end); dsdt_L_theta(2,end)];
dsdt_L_theta = dsdt_L_theta*dti;

dsdt_L = dsdt_L_r + dsdt_L_theta; % numerical left derivative
dsdt_L_plot = [mean(dsdt_L(1,:)); mean(dsdt_L(2,:))];


% numerical derivatives and analytical derivatives all match so far

% for LHS and RHS, the tangent line is in line with the original
% parameterized line. This means that the angle is confirmed to be a
% constant dtheta/2
ds_T = norm(ds_T_dt(:,1)); % numerical norm of top tangent vector
ds_B = norm(ds_B_dt(:,1)); % numerical norm of top tangent vector
ds_R = norm(ds_R_dt(:,1)); % numerical norm of top tangent vector
ds_L = norm(ds_L_dt(:,1)); % numerical norm of top tangent vector

t_hat_T = ds_T_dt/ds_T; % unit tangent top
t_hat_B = ds_B_dt/ds_B; % unit tangent bottom
t_hat_R = ds_R_dt/ds_R;% unit tangent RHS
t_hat_L = ds_L_dt/ds_L; % unit tangent LHS

n_hat_T = [t_hat_T(2,:); -t_hat_T(1,:)]; % outward unit normal top
n_hat_B = [t_hat_B(2,:); -t_hat_B(1,:)]; % outward unit normal bottom
n_hat_R = [t_hat_R(2,:); -t_hat_R(1,:)]; % outward unit normal RHS
n_hat_L = [t_hat_L(2,:); -t_hat_L(1,:)]; % outward unit normal LHS


frac_n_T_dA = n_hat_T*ds_T/dA*dt;
frac_n_B_dA = n_hat_B*ds_B/dA*dt;
frac_n_R_dA = n_hat_R*ds_R/dA*dt;
frac_n_L_dA = n_hat_L*ds_L/dA*dt;

% I have confirmed that the tangent vectors are, in fact, tangent
% I have confirmed that the normal vectors are, in fact, normal


% generate scalar field f:

x_T = (r0+0.5*dr)*cos(theta0+dtheta*t);
y_T = (r0+0.5*dr)*sin(theta0+dtheta*t);

x_B = (r0-0.5*dr)*cos(theta0-dtheta*t);
y_B = (r0-0.5*dr)*sin(theta0-dtheta*t);

x_R = (r0+dr*t)*cos(theta0-0.5*dtheta);
y_R = (r0+dr*t)*sin(theta0-0.5*dtheta);

x_L = (r0-dr*t)*cos(theta0+0.5*dtheta);
y_L = (r0-dr*t)*sin(theta0+0.5*dtheta);

f = @(x,y) 2*x + x.*y + y.^2;
% Evaluate field over parameterized curve:
f_T = f(x_T,y_T);
f_B = f(x_B,y_B);
f_R = f(x_R,y_R);
f_L = f(x_L,y_L);

vect_T = f_T.*frac_n_T_dA;
vect_B = f_B.*frac_n_B_dA;
vect_R = f_R.*frac_n_R_dA;
vect_L = f_L.*frac_n_L_dA;

Ans_T = sum(vect_T,2);
Ans_B = sum(vect_B,2);
Ans_R = sum(vect_R,2);
Ans_L = sum(vect_L,2);

Ans = Ans_T + Ans_B + Ans_R + Ans_L;
dfdx = Ans(1);
dfdy = Ans(2);
disp('Partial f partial x:')
disp(dfdx)
disp('partial f partial y:')
disp(dfdy)




figure
plot(s_R(1,:),s_R(2,:))
hold on
plot(s_L(1,:),s_L(2,:))
hold on
plot(s_T(1,:),s_T(2,:))
hold on
plot(s_B(1,:),s_B(2,:))
hold on
plot([0 5],[0 0],'k','LineWidth',3)
hold on
plot([0 0],[0 5],'k','LineWidth',3)
hold on
% quiver(s_T(1,:),s_T(2,:),t_hat_T(1,:),t_hat_T(2,:))
% hold on
% quiver(s_B(1,:),s_B(2,:),t_hat_B(1,:),t_hat_B(2,:))
% hold on
% quiver(s_R(1,:),s_R(2,:),t_hat_R(1,:),t_hat_R(2,:))
% hold on
% quiver(s_L(1,:),s_L(2,:),t_hat_L(1,:),t_hat_L(2,:))

% quiver(s_T(1,:),s_T(2,:),n_hat_T(1,:),n_hat_T(2,:))
% hold on
% quiver(s_B(1,:),s_B(2,:),n_hat_B(1,:),n_hat_B(2,:))
% hold on
% quiver(s_R(1,:),s_R(2,:),n_hat_R(1,:),n_hat_R(2,:))
% hold on
% quiver(s_L(1,:),s_L(2,:),n_hat_L(1,:),n_hat_L(2,:))
% hold on

quiver(s_T(1,:),s_T(2,:),frac_n_T_dA(1,:),frac_n_T_dA(2,:))
hold on
quiver(s_B(1,:),s_B(2,:),n_hat_B(1,:),n_hat_B(2,:))
hold on
quiver(s_R(1,:),s_R(2,:),n_hat_R(1,:),n_hat_R(2,:))
hold on
quiver(s_L(1,:),s_L(2,:),n_hat_L(1,:),n_hat_L(2,:))
hold on
plot(x_B(2),y_B(2),'r*','MarkerSize',5)

% the direction of each one is correct, so f is being evaluated properly



%hold on
%plot(ds_B_dt(1,:),ds_B_dt(2,:))
%hold on
%plot([0 ds_R_dt(1,:)],[0 ds_R_dt(2,:)])
%hold on
%plot([0 -ds_L_dt(1,:)],[0 -ds_L_dt(2,:)])


% hold on
% plot(n_hat_T(1,:),n_hat_T(2,:),'LineWidth',1.5)
% hold on
% plot(n_hat_B(1,:),n_hat_B(2,:),'LineWidth',1.5)
legend('polar coordinate RHS','polar coordinate LHS','polar coordinate top','polar coordinate bottom','polar tangent top','polar tangent bottom','polar tangent RHS','polar tangent LHS','location','SW')
axis equal







% plot S_R under the polar and then under the cartesian basis to see how
% its look changes


% also compare the parameterization from the new way you were computing the
% parameterization to the old way


% additionally, see if the s_R and s_L are in line with their tangent
% counterparts
