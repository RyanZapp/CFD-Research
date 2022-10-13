clear
close all
clc
%pp
Lx = 2;
Ly = 1;
rho = 1 ;
nu = 0.008926 ;
%t_final = 2 ;
%N = 10000 ;
%J = 250 ;


%mu = rho*nu ;
% Eventually swap h with dx and dy
%h = Ly/(J-1) ;



nx = 81; % must be a perfect square
ny = 81; % must be a perfect square
N = 100;
% Create indexing
imin = 2;
jmin = 2; % Start from 2 to allow for ghost nodes
imax = imin + nx - 1; % I think I need to remove the -1 from this and the line
% below in order to have the correct number of interior points
jmax = jmin + ny - 1;

% Generate mesh
x(imin:imax+1) = linspace(0,Lx,nx+1);
y(jmin:jmax+1) = linspace(0,Ly,ny+1); % Basically, we are doing jmax+1 to add a ghost node
% ans we compansate for that by taking ny+1

% Create mesh sizes
dx = x(imin+1) - x(imin); % Basically just computing the distance between
% two cells and saying that the calue will be replicated for all inter-cell
% distances
dy = y(jmin+1) - y(jmin);
dxi = 1/dx;
dyi = 1/dy;
dxi2 = dxi^2;
dyi2 = dyi^2;
dt = 0.125/nu*dxi^2 ;
% We precompute the division by dx and dy becasue division is a lot more
% computationally expensive than multiplication
%x = linspace(0,L,J) ;
%y = linspace(0,ht,J) ;
%k = -0.5*x + 1;
%t = linspace(0,t_final,N); % There is no point in storing t, only deltaT matters
u = zeros(nx,ny);
v = zeros(nx,ny);
%u(:,imax+1) = neumann
u(:,imin-1)
u(:,imin-1) = y.^2
% doing this essentially initialized the no slip BC for free

% In this region, we need to specify boundary conditions, and then we need
% to specify an initial condition

% Actually, we need to first create the original u vector with the initial
% condition, then we need to append the first and last columns with the
% boundary conditions

% We also need to append the first and last rows with the boundary
% conditions

[u0,v0] = initCond(x,y,dxi,dyi,nx,ny,Lx,Ly);

% u(:,:,1) = u0;
% v(:,:,1) = v0;
% 
% for i = 2:N
%     u(:,1,i) = u0(:,1);
%     v(:,1,i) = v0(:,1);
% end
% 
% 
% % By this line, the I.C. should be det, along with the B.C. for the left,
% % top, and bottom sides
%  % We could just add in another generic loop here to make u keep going
%  % and overriding itself, which would have the effect of the solution
%  % marching forward in time
for j = jmin:jmax
    for i = imin+1:imax % why are we leaving node 2 out?
        v_here = 0.25*(v(i-1,j) + v(i,j) + v(i-1,j) + v(i,j+1));
        us(i,j) = u(i,j) + dt*(nu*((u(i-1,j)-2*u(i,j)+u(i+1,j))*dxi^2 + (u(i,j-1)-2*u(i,j)+u(i,j+1))*dyi^2) ... 
            - (u(i,j)*(u(i+1,j)-u(i-1,j))*0.5*dxi + v_here*(u(i,j+1)-u(i,j-1))*0.5*dyi));
    end
    u(i,imax+1) = u(i,imax); % Neumann BC (might be implmented wrong)
end

for j = jmin+1:jmax
    for i = imin:imax
        u_here = 0.25*(u(i,j-1) + u(i,j) + u(i+1,j-1) + u(i+1.j));
        vs(i,j) = v(i,j) + dt*(nu*((v(i-1,j)-2*v(i,j)+v(i+1,j))*dxi^2 + (v(i,j-1)-v*u(i,j)+v(i,j+1))*dyi^2) ... 
            -(u_here*(v(i+1,j)-v(i-1,j))*0.5*dxi + v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*dyi);
    end
    v(i,jmax+1) = v(i,jmax); % Neumann BC (might be implemented wrong)
end

% for n = 1:N-1
%     for i = 2:J-1
%         for j = 2:J-2
%            u(i,j,n+1) = deltaT*nu/(h^2)*(u(i+1,j,n) + u(i,j+1,n) - 4*u(i,j,n) + u(i-1,j,n) + u(i,j-1,n)) - deltaT*u(i,j,n)/h*(u(i+1,j,n) - u(i,j,n)) - deltaT*v(i,j,n)/h*(u(i,j+1,n) - u(i,j,n)) + u(i,j,n) ;
%            v(i,j,n+1) = deltaT*nu/(h^2)*(v(i+1,j,n) + v(i,j+1,n) - 4*v(i,j,n) + v(i-1,j,n) + v(i,j-1,n)) - deltaT*u(i,j,n)/h*(v(i+1,j,n) - v(i,j,n)) - deltaT*v(i,j,n)/h*(v(i,j+1,n) - v(i,j,n)) + v(i,j,n) ;
% 
%         end
%            u(i,J-1,n+1) = u(i,J-2,n+1); % Impose Neumann B.C. on rhs for x velocity
%            v(i,J-1,n+1) = v(i,J-2,n+1); % Impose Neumann B.C. on rhs for y velocity
%     end
%     
% end

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
