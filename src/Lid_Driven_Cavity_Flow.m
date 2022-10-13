clear
close all
clc
rho = 1 ;
nu = 0.008926 ;
nx = 81;
ny = 81;
Lx = 1;
Ly = 1;
dx = (Lx-0)/nx;
dy = (Ly-0)/ny;
% Predivide and pre-square becasue multiplication is more computationally
% efficient than division
dxi = 1/dx;
dyi = 1/dy;
dxi2 = dxi^2;
dyi2 = dyi^2;

nt = 10; % max number of timesteps

dt = 0.00001*dx*dy/(nu*(dx + dy));% timestep size (determined by CFL condition)
%dt = 0.125*dxi2/nu; % This is a backup timestep size in case the other one
%fails

% Specify location of cell centers (this is where we will be plotting)
xce = ((1:nx)-0.5)*dx; % length(xce) = nx
yce = ((1:ny)-0.5)*dy; % length(yce) = ny

% Create vector that houses spatial x coordinates (starting from bottom
% left corner of physical domain)
xco = (0:nx)*dx;  % length(xco) = nx + 1
yco = (0:ny)*dy; % length(yco) = ny + 1

% initialize matrices for u, v, and p
u = zeros(nx+1,ny+2);
us = zeros(nx+1,ny+2);
v = zeros(nx+2,ny+1);
vs = zeros(nx+2,ny+1);
p = zeros(nx,ny);
R1 = zeros(nx,ny);
R2 = zeros(nx,ny);
uce = zeros(nx,ny);
vce = zeros(nx,ny);
%note that p is the only matrix of the three without ghost nodes

% call the pressure matrix function
% The pressure matrix (i.e. the Laplacian) depends only on the geometry of
% the system, so we can initialize it outside of the loop to increase
% efficiency

L = pressure_matrix(nx,ny,dxi2,dyi2);

% to access the interior points of each matrix, type the following:
% u(2:end-1,2:end-1);
% v(2:end-1,2:end-1);
% p(1:end,1:end);

% at this point, if you have a nonzero initial condition, you would
% populate u and v with that initial condition here

% however, our initial condition is the zeros matrix, so we will proceed.

% Create the meshgrid that you will eventually plot over
[Xce,Yce] = meshgrid(xce,yce);
%figure
[~,h_abs] = contourf(Xce',Yce',sqrt(uce.^2+vce.^2));
xlim([0 Lx]); 
ylim([0 Ly]);

% I need to make the center u and v to make sure that they are plottable
% enter the time iteration:
ubcb = 0;
ubct = 1;
vbcl = 0;
vbcr = 0;

for ii = 1:nt
    % compute B.C.
    u(:,1) = 2*ubcb - u(:,2);         % x-velocity bottom B.C.
    u(:,end) = 2*ubct - u(:,end-1);   % x-velocity top B.C.
    u(1,:) = 0;                       % x-velocity left B.C.
    u(end,:) = 0;                     % x-velocity right B.C.
    v(:,1) = 0;                       % y-velocity bottom B.C.
    v(:,end) = 0;                     % y-velocity top B.C.
    v(1,:) = 2*vbcl - v(2,:);         % y-velocity left B.C.
    v(end,:) = 2*vbcr - v(end-1,:);   % y-velocity right B.C.
    % Repeat for the us and vs matrices
    
    us(:,1) = 2*ubcb - us(:,2);         % x-velocity bottom B.C.
    us(:,end) = 2*ubct - us(:,end-1);   % x-velocity top B.C.
    us(1,:) = 0;                       % x-velocity left B.C.
    us(end,:) = 0;                     % x-velocity right B.C.
    vs(:,1) = 0;                       % y-velocity bottom B.C.
    vs(:,end) = 0;                     % y-velocity top B.C.
    vs(1,:) = 2*vbcl - vs(2,:);         % y-velocity left B.C.
    vs(end,:) = 2*vbcr - vs(end-1,:);   % y-velocity right B.C.
    
    % These matrices are exteremely messed up********
    % note that us is somehow coming out wrong
    % Chris and I need to fix them before proceeding
    
    % solve for u* (u* does not satisfy the continuity equation)
    for j = 2:(ny+1)
        for i = 2:nx
            v_here = 0.25*(v(i+1,j) + v(i+1,j-1) + v(i,j-1) + v(i,j));
            us(i,j) = u(i,j) + dt*(nu*((u(i-1,j)-2*u(i,j)+u(i+1,j))*dxi2 + (u(i,j-1)-2*u(i,j)+u(i,j+1))*dyi2) ... 
            - (u(i,j)*(u(i+1,j)-u(i-1,j))*0.5*dxi + v_here*(u(i,j+1)-u(i,j-1))*0.5*dyi));
        end
        us(:,1) = 2*ubcb - us(:,2);         % x-velocity bottom B.C.
    end
    us(:,end) = 2*ubct - us(:,end-1);   % x-velocity top B.C.
    us(1,:) = 0;                       % x-velocity left B.C.
    us(end,:) = 0;                     % x-velocity right B.C.
    % we put the dirichlet B.C. again becasue the Neumann puts the velocity
    % at the top corners of the box as being equal to the lid velocity,
    % when in reality, we should have a no-slip at the corners.
    % Reintroducing the Dirichlet B.C. fixes this issue
   
    % solve for v* (u* does not satisfy the continuity equation)
    for j = 2:ny
        for i = 2:(nx+1)
            u_here = 0.25*(u(i,j+1) + u(i-1,j) + u(i,j) + u(i-1,j+1));
            vs(i,j) = v(i,j) + dt*(nu*((v(i-1,j)-2*v(i,j)+v(i+1,j))*dxi2 + (v(i,j-1)-2*v(i,j)+v(i,j+1))*dyi2) ... 
            -(u_here*(v(i+1,j)-v(i-1,j))*0.5*dxi + v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*dyi));
        end
        vs(1,:) = 2*vbcl - vs(2,:);         % y-velocity left B.C.
        vs(end,:) = 2*vbcr - vs(end-1,:);   % y-velocity right B.C.
    end
    vs(:,1) = 0;                       % y-velocity bottom B.C.
    vs(:,end) = 0;                     % y-velocity top B.C.
    
    for j = 2:(ny+1)
        for i = 1:nx
            R1(i,j-1) = (us(i+1,j) - us(i,j))*0.5*dxi;
        end
    end
    
    for j = 1:ny
        for i = 2:(nx+1)
            R2(i-1,j) = (vs(i,j+1) - vs(i,j))*0.5*dyi;
        end
    end
    R = -rho/dt*(R1 + R2);
    R = R(:);
    L(1,:) = 0;
    L(1,1) = 1; % for some reason, this stopped the matrix from being singluar....now, the only issue is that, for some reason
    % the velocity for u is degenerating at the top boundary to zero, so
    % the cavity flow quickly comes to a halt
    % I do not know if the entries of R will match up with the
    % corresponding entries of the laplacian
    p = L\R;
    p = reshape(p,nx,ny); 
    LL = inv(L);
    
    % after I fix L, I just need to add the pressure derivatives back into
    % the solution
    
    for j = 2:(ny+1)
        for i = 2:nx
            u(i,j) = us(i,j) - 1/rho*(p(i,j-1) - p(i-1,j-1))*0.5*dxi; % I think this one is correct, but it would not hurt to double check
        end
       
    end
    u(:,1) = 2*ubcb - u(:,2);
    u(:,end) = 2*ubct - u(:,end-1);
 
    u(1,:) = 0;
    u(end,:) = 0;
    for j = 2:ny
        for i = 2:(nx+1)
            v(i,j) = vs(i,j) - 1/rho*(p(i-1,j) - p(i-1,j-1))*0.5*dyi; % I think this is right, but I should double check it
        end
    end
    % something weird is happening in the middle row of the v matrix
    for j = 1:ny
        for i = 1:nx
            uce(i,j) = 0.5*(u(i,j+1) + u(i+1,j+1));
            vce(i,j) = 0.5*(v(i+1,j+1) + v(i+1,j));
        end
        v(1,:) = 2*vbcl - v(2,:);         % y-velocity left B.C.
        v(end,:) = 2*vbcr - v(end-1,:);   % y-velocity right B.C.
    end
    v(:,1) = 0;                     
    v(:,end) = 0;                    
  
    h_abs.ZData = sqrt(uce.^2+vce.^2);
    drawnow
% There are more efficient ways to do this code, but let us just make sure
% it works for now

% The pressure matrix has issues
% I am not sure why for the v matrix, we initialized the Neumann
% boundary conditions in the loop...but in the u matrix, we did the bottom
% Neumann B.C. in the loop and the top Neumann B.C. out of the loop

% I honestly think that both can go outside of the loop
end