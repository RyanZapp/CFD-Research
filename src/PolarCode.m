%% Problem 1
clear
close all
clc
% Shrinkable parameters
N = 10000; % makes parameterization more accurate
dr = 0.0000001; % shrinking this is analogous to taking the limit as dr -> 0
dtheta = 0.0000001; % shrinking this is analogous to taking the limit as dtheta -> 0

b = 0.5;
a = -0.5;
t = linspace(a,b,N)';
dt = (b-a)/(N-1); % we can compute dt rather than using the diff command on t
dti = 1/dt;

theta0 = pi/4;
r0 = 5;

dA = r0*dr*dtheta;


r_hat = [cos(theta0); sin(theta0)];
theta_hat = [-sin(theta0); cos(theta0)];



s_T_r = (r0+0.5*dr)*cos(dtheta*t);
s_T_theta = (r0+0.5*dr)*sin(dtheta*t);
%s_T = s_T_r + s_T_theta;

s_B_r = (r0-0.5*dr)*cos(dtheta*t);
s_B_theta = -(r0-0.5*dr)*sin(dtheta*t);
%s_B = s_B_r + s_B_theta;

s_R_r = (r0+dr*t)*cos(dtheta/2);
s_R_theta = -(r0+dr*t)*sin(dtheta/2);
%s_R = s_R_r + s_R_theta;


s_L_r = (r0-dr*t)*cos(dtheta/2);
s_L_theta = (r0-dr*t)*sin(dtheta/2);
%s_L = s_L_r + s_L_theta;
% first row is the x component and 2nd row is the y component

% my parameterizations look quite good

% Compute analytical derivatives:
ds_T_rdt = -dtheta*(r0+0.5*dr)*sin(dtheta*t);
ds_T_thetadt = dtheta*(r0+0.5*dr)*cos(dtheta*t);
%ds_T_dt = ds_T_rdt + ds_T_thetadt;

ds_B_rdt = -dtheta*(r0-0.5*dr)*sin(dtheta*t);
ds_B_thetadt = -dtheta*(r0-0.5*dr)*cos(dtheta*t);
%ds_B_dt = ds_B_rdt + ds_B_thetadt;

ds_R_rdt = dr*cos(0.5*dtheta);
ds_R_thetadt = -dr*sin(0.5*dtheta);
%ds_R_dt = ds_R_rdt + ds_R_thetadt;
for i = 1:N
    ds_R_rdt(i) = ds_R_rdt(1);
    ds_R_thetadt(i) = ds_R_thetadt(1);
end
ds_R_rdt = ds_R_rdt';
ds_R_thetadt = ds_R_thetadt';

ds_L_rdt = -dr*cos(0.5*dtheta);
ds_L_thetadt = -dr*sin(0.5*dtheta);
%ds_L_dt = ds_L_rdt + ds_L_thetadt;
for i = 1:N
    ds_L_rdt(i) = ds_L_rdt(1);
    ds_L_thetadt(i) = ds_L_thetadt(1);
end
ds_L_rdt = ds_L_rdt';
ds_L_thetadt = ds_L_thetadt';


% numerical derivatives and analytical derivatives all match so far

% for LHS and RHS, the tangent line is in line with the original
% parameterized line. This means that the angle is confirmed to be a
% constant dtheta/2
ds_T = dtheta*(r0+0.5*dr); % theoretical norm of top tangent vector
ds_B = dtheta*(r0-0.5*dr); % theoretical norm of top tangent vector
ds_R = dr; % theoretical norm of top tangent vector
ds_L = dr; % theoretical norm of top tangent vector

t_hat_T = 1/ds_T*[ds_T_rdt,ds_T_thetadt]; % unit tangent top
t_hat_B = 1/ds_B*[ds_B_rdt,ds_B_thetadt]; % unit tangent bottom
t_hat_R = 1/ds_R*[ds_R_rdt,ds_R_thetadt];% unit tangent RHS
t_hat_L = 1/ds_L*[ds_L_rdt,ds_L_thetadt]; % unit tangent LHS

n_hat_T = [t_hat_T(:,2),-t_hat_T(:,1)]; % outward unit normal top
n_hat_B = [t_hat_B(:,2),-t_hat_B(:,1)]; % outward unit normal bottom
n_hat_R = [t_hat_R(:,2),-t_hat_R(:,1)]; % outward unit normal RHS
n_hat_L = [t_hat_L(:,2),-t_hat_L(:,1)]; % outward unit normal LHS


frac_n_T_dA = n_hat_T*ds_T*dt/dA;
frac_n_B_dA = n_hat_B*ds_B*dt/dA;
frac_n_R_dA = n_hat_R*ds_R*dt/dA;
frac_n_L_dA = n_hat_L*ds_L*dt/dA;

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

Ans_T = sum(vect_T);
Ans_B = sum(vect_B);
Ans_R = sum(vect_R);
Ans_L = sum(vect_L);

Term1 = Ans_T(1) + Ans_B(1);
Term_question = Ans_R(1) + Ans_L(1);
Term_null = Ans_T(2) + Ans_B(2);
Term2 = Ans_R(2) + Ans_L(2);
disp('1/r*df/dtheta:')
disp(Term2)
disp('Should be a zero sum (arises from theta_hat coefficient in the top and bottom integral):')
disp(Term_null)
disp('Coefficient corresponding to r_hat in the left and right integral sum:')
disp(Term_question)
disp('If we subtract 6.4149 from this term it should be equal to df/dr:')
disp(Term1)
disp('df/dr:')
Term11 = Term1 + Term_question;
disp(Term11)
Term22 = -2*sin(theta0) - r0*sin(theta0)*sin(theta0) + r0*cos(theta0)*cos(theta0) + 2*r0*sin(theta0)*cos(theta0)

Term111 = 2*cos(theta0) + 2*r0*cos(theta0)*sin(theta0) + 2*r0*sin(theta0)*sin(theta0)
% When I convert to cartesian coordinates, I get the gradient vector at the
% x,y coordinate point corresponding to r0,theta0
% Ans = Ans_T + Ans_B + Ans_R + Ans_L;
% dfdx = Ans(1);
% dfdy = Ans(2);
% disp('Partial f partial x:')
% disp(dfdx)
% disp('partial f partial y:')
% disp(dfdy)




% figure
% plot(s_R(1,:),s_R(2,:))
% hold on
% plot(s_L(1,:),s_L(2,:))
% hold on
% plot(s_T(1,:),s_T(2,:))
% hold on
% plot(s_B(1,:),s_B(2,:))
% hold on
% plot([0 5],[0 0],'k','LineWidth',3)
% hold on
% plot([0 0],[0 5],'k','LineWidth',3)
% hold on
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

% quiver(s_T(1,:),s_T(2,:),frac_n_T_dA(1,:),frac_n_T_dA(2,:))
% hold on
% quiver(s_B(1,:),s_B(2,:),n_hat_B(1,:),n_hat_B(2,:))
% hold on
% quiver(s_R(1,:),s_R(2,:),n_hat_R(1,:),n_hat_R(2,:))
% hold on
% quiver(s_L(1,:),s_L(2,:),n_hat_L(1,:),n_hat_L(2,:))
% hold on
% plot(x_B(2),y_B(2),'r*','MarkerSize',5)

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
% legend('polar coordinate RHS','polar coordinate LHS','polar coordinate top','polar coordinate bottom','polar tangent top','polar tangent bottom','polar tangent RHS','polar tangent LHS','location','SW')
% axis equal

% a nice thing is that the previous code and analytical computations
% confirmed that my parameterization, tangent vectors, and normal vectors
% are all correct...this means that there must be something I am messing up
% with the algebra





% plot S_R under the polar and then under the cartesian basis to see how
% its look changes


% also compare the parameterization from the new way you were computing the
% parameterization to the old way


% additionally, see if the s_R and s_L are in line with their tangent
% counterparts

%% Problem 2
clear
close all
clc
% Initialize constant values
R = 2;
r = 1;
N = 20000;
a = 0;
b = 2*pi;
[theta,phi] = meshgrid(linspace(a,b,N));
% Declare parameterization
x = R*cos(theta) + r*cos(theta).*cos(phi);
y = R*sin(theta) + r*sin(theta).*cos(phi);
z = r*sin(phi);
% Compute partial derivatives
x_theta = -R*sin(theta) - r*sin(theta).*cos(phi);
x_phi = -r*cos(theta).*sin(phi);
y_theta = R*cos(theta) + r*cos(theta).*cos(phi);
y_phi = -r*sin(theta).*sin(phi);
J_det = (x_theta.*y_phi - x_phi.*y_theta);

h = z.^3;
dTheta = (b-a)/(N-1);
dPhi = dTheta;
dA = dTheta*dPhi;
% Compute flux
flux_vect = h.*J_det*dA;
flux = sum(flux_vect(:));
disp('Numerical Flux:')
disp(flux)

% Below is a method that strictly uses differential forms, but the answer
% appears to be half of what the correct answer is...this leads me to
% believe that the method below is viable, but is flawed somehow

% Compute differentials numerically
dx = diff(x,1,2);
dy = diff(y,1,1);
% Fix sizing of matrixes
dx(:,end+1) = dx(:,end);
dy(end+1,:) = dy(end,:);
% Compute wedge product
dx_wedge_dy = 2*dx.*dy;
% Compute flux
Flux_vect = h.*dx_wedge_dy;


Flux = sum(Flux_vect(:));
disp('Numerical Flux using the language of forms:')
disp(Flux)


