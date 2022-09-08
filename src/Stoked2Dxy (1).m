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

h = ht/(J-1) ;
deltaT = 0.125/nu*h^2 ;

x = linspace(0,L,J) ;
y = linspace(0,ht,J) ;
k = -0.5*x + 1;
t = linspace(0,t_final,N);

u = zeros(J,J,N);
v = zeros(J,J,N);

%u0 = zeros(J,J);

u(2:(J-1),1,:) = 0.5*ones((J-2),1,N) ; %if error, start trouble shooting here. need to be able to enforce condition for all N
%u(2:(J-1),J,:) = ones((J-2),1,N) ; %rationalized using volumetric flow for uniform cross sectional area Qin = Qout
%direchlet/periodic condition
M = 0.5*ones(J,J);
for i = 2:J-1
    for j = 2:J-2
    M(i,j) = k(j);
    end
    M(i,J-1) = M(i,J-2);
end
M(1,:) = 0;
M(J,:) = 0;
M(:,J) = 0;
u(:,:,1) = M;
v(:,:,1) = IC;
v(:,1,:) = % we want to set velocity I.C. as a parabola that doesnt have much magnitude
% The b.c. should be the first time step of that parabola...not sure where
% the file I made for it went tho
% set neumann
% The solutions exist and are smooth for I.C. = 0, I.C. = polynomial, I.C.
% = sine wave that has reverse and forward flow, and much more...I.e. this
% solution is very stable (letting u0 = 0 is akin to having a system at
% rest and then introducing a constant pressure gradient.) That being said,
% it implies that we can just use 0 I.C. if we need to in our 2D testerd

for n = 1:N-1
    for i = 2:J-1
        for j = 2:J-2
           u(i,j,n+1) = deltaT*nu/(h^2)*(u(i+1,j,n) + u(i,j+1,n) - 4*u(i,j,n) + u(i-1,j,n) + u(i,j-1,n)) - deltaT*u(i,j,n)/h*(u(i+1,j,n) - u(i,j,n)) - deltaT*v(i,j,n)/h*(u(i,j+1,n) - u(i,j,n)) + u(i,j,n) ;
           v(i,j,n+1) = deltaT*nu/(h^2)*(v(i+1,j,n) + v(i,j+1,n) - 4*v(i,j,n) + v(i-1,j,n) + v(i,j-1,n)) - deltaT*u(i,j,n)/h*(v(i+1,j,n) - v(i,j,n)) - deltaT*v(i,j,n)/h*(v(i,j+1,n) - v(i,j,n)) + v(i,j,n) ;
           %u(j,i,n+1) = deltaT*(deltaP/(rho*L) + nu/(h^2)*(u(j+1,i,n) + u(j,i+1,n) - 4*u(j,i,n) + u(j-1,i,n) + u(j,i-1,n)) - u(j,i,n)/h*(u(j+1,i,n) - u(j,i,n)) ) + u(j,i,n) ;

        end
           u(i,J-1,n+1) = u(i,J-2,n+1);
           v(i,J-1,n+1) = v(i,J-2,n+1);
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
