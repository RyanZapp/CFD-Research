clear
close all
clc
x = linspace(-2*pi,2*pi);
y = linspace(0,4*pi);
[X,Y] = meshgrid(x,y);
Z = sin(X) + cos(Y);
Z(:,1) = -5; % this should be the bottom after the transform
Z(:,2) = -5;
Z(:,3) = -5;
Z(:,4) = -5;
Z(:,end) = 1; % this should be top after the transform
Z(1,:) = 5; % this should be left after the transform
Z(2,:) = 5;
Z(3,:) = 5;
Z(4,:) = 5;
Z(end,:) = 0; % This should be right after the transform
Z(end-1,:) = 0;
Z(end-2,:) = 0;
Z(end-3,:) = 0;
Z(end-4,:) = 0;
KX = X';
KY = Y';
figure
subplot(1,2,1)
%contourf(X,Y,Z,10)
%s1 = pcolor(X,Y,Z)
%set(s1, 'EdgeColor', 'none');
imagesc(Z)
%set(gca,'YDir','normal')
axis equal
title("Imagesc Before transform")
subplot(1,2,2)
imagesc(Z')
set(gca,'YDir','normal') % doing this line flips the y direction properly to give us the proper output image
%contourf(X,Y,Z,10)
%s2 = pcolor(X',Y',Z)
%set(s2, 'EdgeColor', 'none');
axis equal
title("Imagesc After transform")

figure
subplot(1,2,1)
contourf(X,Y,Z,10)
axis equal
title("Contourf Before transform")
subplot(1,2,2)
contourf(X',Y',Z,10)
axis equal
title("Contourf After transform")



% After careful analysis, it is apparent that contourf, pcolor, and imagesc
% are all valuable methods of plotting the outputs and can all be made to
% work