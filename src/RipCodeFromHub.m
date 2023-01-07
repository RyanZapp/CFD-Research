clear
close all
clc

Re = 300; % Reynolds number
nt = 300; % max time steps (お試しで少しだけ)
Lx = 1; Ly = 1; % domain size
Nx = 51; Ny = 50; % Number of grids
dt = 0.01; % time step;

% Grid size (Equispaced)
dx = Lx/Nx;
dy = Ly/Ny;
% Coordinate of each grid (cell center)
xce = ((1:Nx)-0.5)*dx;
yce = ((1:Ny)-0.5)*dy;
% Coordinate of each grid (cell corner)
xco = (0:Nx)*dx;
yco = (0:Ny)*dy;

u = zeros(Nx+1,Ny+2); % velocity in x direction (u)
v = zeros(Nx+2,Ny+1); % velocity in y direction (v)
p = zeros(Nx,Ny); % pressure (lagurange multiplier)

uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2; % u at cell center % the first term goes u(1,2), u(1,3), u(1,4), u(1,5) and then it gets to u(1,end-1)
% then it goes u(2,2), u(2,3), u(2,4), and so on, until it gets to
% u(2,end-1). It keeps incrementing this way until it gets to
% u(end-1,end-1)... the 2nd term works this way as well..
% it is essentially the same as doing the following loop
% uce = zeros(Nx,Ny);
% for i = 1:(end-1)
    %for j = 2:(end-1)
        % uce(i,j-1) = (u(i,j)+u(i+1,j))*0.5;
    %end
%end
vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2; % v at cell center
% The same idea for uce generation applies to vce
[Xce,Yce] = meshgrid(xce,yce); % cell centerの座標グリッド
[~,h_abs] = contourf(Xce',Yce',sqrt(uce.^2+vce.^2));

xlim([0 Lx]); ylim([0 Ly]);

for ii = 1:nt
    bctop = 1; % Top lid u
    u(:,1) = -u(:,2); v(:,1) = 0;             %bottom
    u(:,end) = 2*bctop-u(:,end-1);  v(:,end) = 0;  %top
    u(1,:) = 0;    v(1,:) = -v(2,:);             %left
    u(end,:) = 0;  v(end,:) = -v(end-1,:);    %right
    % They put dirichlet B.C. after Neumann (like we did) in order to make
    % sure that the values at the corners of the top remain as zero
    Lux = (u(1:end-2,2:end-1)-2*u(2:end-1,2:end-1)+u(3:end,2:end-1))/dx^2; % nx-1 * ny
    Luy = (u(2:end-1,1:end-2)-2*u(2:end-1,2:end-1)+u(2:end-1,3:end))/dy^2; % nx-1 * ny
    Lvx = (v(1:end-2,2:end-1)-2*v(2:end-1,2:end-1)+v(3:end,2:end-1))/dx^2; % nx * ny-1
    Lvy = (v(2:end-1,1:end-2)-2*v(2:end-1,2:end-1)+v(2:end-1,3:end))/dy^2; % nx * ny-1
    % shown above are the 2nd derivatives. I have not yet done an in depth
    % check of if they match up with what we have
    
    uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2;
    uco = (u(:,1:end-1)+u(:,2:end))/2; % what are the co values?
    vco = (v(1:end-1,:)+v(2:end,:))/2;
    vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2;
    % I think the lines above just reinstitute the uce and vce, but this
    % time, they carry over the boundary conditions from u and v
    % 2. multiply
    uuce = uce.*uce;
    uvco = uco.*vco;
    vvce = vce.*vce;
    
    % 3-1. get derivative for u
    Nu = (uuce(2:end,:) - uuce(1:end-1,:))/dx;
    Nu = Nu + (uvco(2:end-1,2:end) - uvco(2:end-1,1:end-1))/dy;
    % 3-2. get derivative for v
    Nv = (vvce(:,2:end) - vvce(:,1:end-1))/dy;
    Nv = Nv + (uvco(2:end,2:end-1) - uvco(1:end-1,2:end-1))/dx;
    
    % Get intermidiate velocity
    u(2:end-1,2:end-1) = u(2:end-1,2:end-1) + dt*(-Nu + (Lux+Luy)/Re);
    v(2:end-1,2:end-1) = v(2:end-1,2:end-1) + dt*(-Nv + (Lvx+Lvy)/Re); % solving for us and vs, but they just call them u and v here
    
    % velocity correction
    % RHS of pressure Poisson eq.
    b = ((u(2:end,2:end-1)-u(1:end-1,2:end-1))/dx ...
        + (v(2:end-1,2:end)-v(2:end-1,1:end-1))/dy); % this is just taking the sum of partial derivatives as the divergence
    % why dont they scale b by rho and dt???
    
    % Solve for p (using cosine transform, faster)
    p = solvePoissonEquation_2dDCT(b,Nx,Ny,dx,dy);
    
    % Direct method
    % p = solvePoissonEquation_direct(b,Nx,Ny,dx,dy);
    
    % The new divergent free velocity field
    u(2:end-1,2:end-1) = u(2:end-1,2:end-1) -  (p(2:end,:)-p(1:end-1,:))/dx;
    v(2:end-1,2:end-1) = v(2:end-1,2:end-1) -  (p(:,2:end)-p(:,1:end-1))/dy;
    
    
    
     % get velocity at the cell center (for visualization)
    uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2;
    vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2;
    
    % update plot
    h_abs.ZData = sqrt(uce.^2+vce.^2);
    drawnow
end
    
    
