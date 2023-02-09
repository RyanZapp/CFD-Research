% Define the geometry of the domain and the cylinder
L = 1; % length of domain
H = 1; % height of domain
r = 0.25; % radius of cylinder
xCylinder = 0.5; % x-coordinate of center of cylinder
yCylinder = 0.5; % y-coordinate of center of cylinder

% Discretize the domain into a grid of points
nx = 50; % number of grid points in the x-direction
ny = 50; % number of grid points in the y-direction
dx = L/(nx-1); % grid spacing in the x-direction
dy = H/(ny-1); % grid spacing in the y-direction
[X,Y] = meshgrid(0:dx:L,0:dy:H);

% Define the boundary conditions for the problem
% Set the velocity at the inlet to 1 in the x-direction
uInlet = ones(ny,1);
vInlet = zeros(ny,1);

% Set the velocity at the outlet to 0 in the x-direction
uOutlet = zeros(ny,1);
vOutlet = zeros(ny,1);

% Set the no-slip condition on the walls and cylinder
uWalls = zeros(nx,1);
vWalls = zeros(nx,1);
uCylinder = zeros(size(X));
vCylinder = zeros(size(X));
u = zeros(nx,ny);
v = zeros(nx,ny);
% Iterate over the grid points and use the finite difference method to 
% approximate the derivatives
for i = 2:nx-1
    for j = 2:ny-1
        % Calculate the x- and y-velocities at the current point
        u(i,j) = (u(i+1,j)-u(i-1,j))/(2*dx); % use finite difference method to approximate du/dx
        v(i,j) = (v(i,j+1)-v(i,j-1))/(2*dy); % use finite difference method to approximate dv/dy
    end
end

% Use an iterative method to solve for the velocity and pressure fields
maxIter = 1000; % maximum number of iterations
tol = 1e-6; % tolerance for convergence
for k = 1:maxIter
    % Update the velocity field using the Jacobi method
    for i = 2:nx-1
        for j = 2:ny-1
            u(i,j) = (u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1))/4;
            v(i,j) = (v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1))/4;
        end
    end
    
    % Check for convergence
    if norm(u-uOld)/norm(uOld) < tol && norm(v-vOld)/norm(vOld) < tol
        break;
    end
    
    % Update the old velocity field
    uOld = u;
    vOld = v;
end

% Visualize the results
quiver(X,Y,u,v)
contourf(X,Y,p)
