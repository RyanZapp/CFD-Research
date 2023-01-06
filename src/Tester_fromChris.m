%% Numerical Gradient Calucation
clear
close all
clc
%For a Quadrilateral
% deltax = 0.000001 ; 
% deltay = 0.000001 ;

r = 0.000001 ;

N = 200000 ;
t = linspace(0,2*pi,N+1) ;

%Parameterization

%Point to Evaluate The Gradient At
x0 = 1 ; 
y0 = 5 ; 

x = x0 + r*cos(t) ;
y = y0 + r*sin(t) ;

a = 3 ;
b = 1 ;
c = 1 ;
d = 2 ;

A = [((a*d - b*c)*pi*r^2) (pi*r^2)] ; 
k = 1 ;
% k = 1 : Area with jacobian mapping
% k = 2 : Area without jacobian mapping

%New coordinate system
u = a.*x + b.*y  ;
v = c.*x + d.*y  ;

%Evaluation of point in uv
u0 = a.*x0 + b.*y0 ;
v0 = c.*x0 + d.*y0 ;

%Function
f = @(u,v) u.*v ;
f1 = @(x,y) (a*x + b*y).*(c*x + d*y);
dx = diff(x);
dx(end+1) = dx(end);
dy = diff(y);
dy(end+1) = dy(end);
ds1 = sqrt(dx.^2 + dy.^2);
nx = dy./ds1;
ny = -dx./ds1;
delx = sum((f1(x,y)/A(2)).*nx.*ds1);
dely = sum((f1(x,y)/A(2)).*ny.*ds1);
detP = a*d-b*c;
eu = 1/detP*[d;-c];
ev = 1/detP*[-b;a];
%Differentials and Arc Length
du = diff(u) ;
dv = diff(v) ; 
ds = sqrt(du.^2 + dv.^2) ;

t(end) = [] ;
u(end) = [] ;
v(end) = [] ;

%Components of the Normal Vector (Clockwise 90 deg rotation of the tangent
%vector)
nu = dv./ds ;
nv = -du./ds ;

delu = sum((f(u,v)/A(k)).*nu.*ds) ;
delv = sum((f(u,v)/A(k)).*nv.*ds) ;

disp('Numerical Gradient')
disp("Numeric f_u")
disp(delu)
disp("Numeric f_v")
disp(delv)
e1 = [1;0];
e2 = [0;1];
%Analytical Solution
adelu = v0 ;
adelx = 2*a*c*x0 + (a*d+b*c)*y0;
adely = (a*d+b*c)*x0 + 2*b*d*y0;
adelv = u0 ;
agraduv = adelu*eu + adelv*ev;
deluv = delu*eu + delv*ev;
agradxy = adelx*e1 + adely*e2;
delxy = delx*e1 + dely*e2;
gradientPx = (a*d - b*c)*(d/(d^2+c^2)*adelu - d/(a*c+d*b)*adelv + b/(a*c+d*b)*adelu - b/(a^2+b^2)*adelv);
gradientPx = (a*d-b*c)*(d/(d^2+c^2)*adelu - b/(a^2+b^2)*adelv);
disp('Analytical Gradient')
disp("Analytical f_u")
disp(adelu)
disp("Analytical f_v")
disp(adelv) 


disp('Cartesian Analytical Value:')
disp(agradxy)
disp('Cartesian Numerical Value:')
disp(delxy)
disp('UV Analytical Value:')
disp(agraduv)
disp('UV Numerical Value:')
disp(deluv)
disp('Alternative for gradient:')
disp(gradientPx)