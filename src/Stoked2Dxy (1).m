clear
close all
clc
%pp
L = 2;
ht = 1;
rho = 1 ;

t_final = 2 ;
N = 10000 ;
J = 250 ;

nu = 0.008926 ;
mu = rho*nu ;
% Eventually swap h with dx and dy
h = ht/(J-1) ;
deltaT = 0.125/nu*h^2 ;

x = linspace(0,L,J) ;
y = linspace(0,ht,J) ;
k = -0.5*x + 1;
t = linspace(0,t_final,N);

u = zeros(J,J,N);
v = zeros(J,J,N);

[u0,v0] = initCond(x,y,h,J,L,ht);

u(:,:,1) = u0;
v(:,:,1) = v0;

for i = 2:N
    u(:,1,i) = u0(:,1);
    v(:,1,i) = v0(:,1);
end
% By this line, the I.C. should be det, along with the B.C. for the left,
% top, and bottom sides

for n = 1:N-1
    for i = 2:J-1
        for j = 2:J-2
           u(i,j,n+1) = deltaT*nu/(h^2)*(u(i+1,j,n) + u(i,j+1,n) - 4*u(i,j,n) + u(i-1,j,n) + u(i,j-1,n)) - deltaT*u(i,j,n)/h*(u(i+1,j,n) - u(i,j,n)) - deltaT*v(i,j,n)/h*(u(i,j+1,n) - u(i,j,n)) + u(i,j,n) ;
           v(i,j,n+1) = deltaT*nu/(h^2)*(v(i+1,j,n) + v(i,j+1,n) - 4*v(i,j,n) + v(i-1,j,n) + v(i,j-1,n)) - deltaT*u(i,j,n)/h*(v(i+1,j,n) - v(i,j,n)) - deltaT*v(i,j,n)/h*(v(i,j+1,n) - v(i,j,n)) + v(i,j,n) ;

        end
           u(i,J-1,n+1) = u(i,J-2,n+1); % Impose Neumann B.C. on rhs for x velocity
           v(i,J-1,n+1) = v(i,J-2,n+1); % Impose Neumann B.C. on rhs for y velocity
    end
    
end

figure
%[x,y] = meshgrid(linspace(0,L,J),linspace(0,ht,J));
%for i = 1:100:5000
s = pcolor(u(:,:,10000)) ;
set(s,'EdgeColor','none') ;
s.FaceColor ='interp' ; 
    %g = surf(x,y,u(:,:,1000),'EdgeAlpha',0.2,'FaceAlpha',1.0);
    colorbar
    %set(g,'LineStyle','none')
    ylabel('Height of channel')
    xlabel('Length of channel')
    title('2D XY Channel flow')
   % hold on
    %plot(vx,y1,'ro')
    %legend('Numerical Transient Solution','Analytical SS solution','location','NE')
   % hold off
    %grid on
    %axis([0 8 0 h]) % Change this when you change your IC
  %  pause(0.00000001)
%end
%[y,z] = meshgrid(linspace(0,ht,J),linspace(0,ht,J));
% We have our surface function u
% g = surf(y,z,u(:,:,1),'EdgeAlpha',0.2,'FaceAlpha',1.0);
% set(g,'LineStyle','none')
% colorbar
% ylabel('Height of channel')
% xlabel('Width of channel')
% title('1D transient Laminar flow with edge affects accounted for')
