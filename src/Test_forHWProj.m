clear
close all
clc
% shrinkable parameters
N = 1000;
dx = 1;
dy = 1;
t = linspace(-1,1,N);
h = 0.001;

x0 = 0;
y0 = 0;

a = 3;
b = 1;
c = 1;
d = 2;

xT = h;
yT = h*t;
\]
,jn b
dA = dx*dy;
det_P = a*d-b*c;

eu = 1/det_P*[d; -c];
ev = 1/det_P*[-b; a];
u0 = a*x0 + b*y0;
v0 = c*x0 + d*y0;

k = ones(1,N);
s_R = [h*k;h*t];
uT = a*s_R(1,:) + b*s_R(2,:);
vT = c*s_R(1,:) + d*s_R(2,:);
duT = diff(uT);
duT(end+1) = duT(end);
dvT = diff(vT);
dvT(end+1) = dvT(end);
s_Ru = a*h + b*h*t;
r_Rv = c*h + d*h*t;
ds1 = sqrt(s_R(1,:).^2 + s_R(2,:).^2);
sum(ds1)
ds2 = 1/det_P*sqrt((c^2 + d^2).*duT.^2 + (a^2 + b^2).*dvT.^2 - 2*(a*c + d*b)*duT.*dvT);
sum(ds2)