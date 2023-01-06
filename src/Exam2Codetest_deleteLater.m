clear
close all
clc
N = 10;
% Cartesian Pre-Processing
[x,y] = meshgrid(linspace(-2,2,N),linspace(-2,2,N));
[xx,yy] = meshgrid(linspace(-2,2,10),linspace(-2,2,10));

h1 = (2 + 2)/(N-1);
f_cartesian = 2*x + x.*y + y.^2;
[fx,fy] = gradient(f_cartesian,h1);
fx_interp = interp2(x,y,fx,xx,yy,'cubic');
fy_interp = interp2(x,y,fy,xx,yy,'cubic');

% Polar Pre-Processing
[r,theta] = meshgrid(linspace(0.001,2,N),linspace(0,2*pi,N));
[rr,tt] = meshgrid(linspace(0.001,2,10),linspace(0,2*pi,10));

er = [cos(theta); sin(theta)];
etheta = [-sin(theta); cos(theta)];

h2_r = (2-0.001)/(N-1);
h2_theta = (2*pi-0)/(N-1);
f_polar = 2*r.*cos(theta) + r.^2.*cos(theta).*sin(theta) + r.^2.*sin(theta).^2;
[fr,ftheta] = gradient(f_polar,h2_r,h2_theta);
FRR = 2*cos(theta) + 2*r.*cos(theta).*sin(theta) + 2*r.*sin(theta).^2;
ftheta_overR = 1./r.*ftheta;
grad_f_polar_x = fr*er(1) + ftheta_overR*etheta(1);
grad_f_polar_y = fr*er(2) + ftheta_overR*etheta(2);
fr_interp = interp2(r,theta,grad_f_polar_x,rr,tt,'cubic');
ftheta_interp = interp2(r,theta,grad_f_polar_y,rr,tt,'cubic');

% Post-Processing
figure
quiver(xx,yy,fx_interp,fy_interp,0.8,'LineWidth',1)
title('Gradient of f in cartesian coordinates')
figure
quiver(rr,tt,fr_interp,ftheta_interp,0.8,'LineWidth',1)
title('Gradient of f in polar coordinates')
