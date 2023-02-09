% Numerical CFD Solver
clear
close all
clc

%% ALGRORITHM OVERVIEW
% Parameters -> Indexing -> Intialization -> Pressure Matrix Generatrion ->
% Laplacian Operator (and its inverse) -> Time Loop {Apply B.C. -> Predictor
% Step -> RHS of Pressure Poisson EQ -> Solve Poisson -> Corrector Step I
% -> Corrector Step II -> Plot Frame for velocity and pressure field ->} END
%(P. 11 OWKES)

% Outline and inspiration from "A Guide To Writing Your First CFD Solver"
% By Mark Owkes, 2017 (See RESOURCES in GITHUB for all EQ and PAGES referenced)

%% Parameters 

% Fluid Intrinsic Properties (This sets units for problem)
rho = 1 ; % currently in (g/cm^3)
nu = 0.008926 ; % currently in (cm^2/s)
%nu = 1/300;
% Next up, I need to test with Lx=1.5 and nu=1/300 with exact matches to
% rip code from hub in order to accurately gage the match between my code
% and theirs
% Calculation mode for Ru and Rv (Default is F = 1)
F = 2;
H = 2;
% R matches b exactly when F = 2 

% Initial Boundary Conditions, B.C. (u and v direction)
ubcb = 0;
ubct = 1;
vbcl = 0;
vbcr = 0;

%% "Container" Space Discretization (xy-plane cross-section) - Indexing

% Number of cells
n = 100; 

% Number of cells in x and y axis
nx = n+20;
ny = n;

% Dimensions of spatial domain
Lx = 1.5;
Ly = 1;

% "step" size in spatial domain
dx = (Lx-0)/nx;
dy = (Ly-0)/ny;

% pre-compute relevant quotients and products
dxi = 1/dx;
dyi = 1/dy;
dxi2 = dxi^2;
dyi2 = dyi^2;

%% Time Domain Discretization - Indexing 

% Max time steps
nt = 500; 

dt = 0.001;
dti = 1/dt;

%% Cell Center Location - Indexing

xce = ((1:nx)-0.5)*dx; % length(xce) = nx
yce = ((1:ny)-0.5)*dy; % length(yce) = ny

%% Spatial Coordinate Vectors (from bottom left corner of domain) - Indexing

xco = (0:nx)*dx;  % length(xco) = nx + 1
yco = (0:ny)*dy; % length(yco) = ny + 1

delta = 0.05;
[xx,yy] = meshgrid(0:delta:Lx,0:delta:Ly); 

%% Velocity (u,v) and Pressure (p) Matrix Initialization 

% u (x-component velocity)
u = zeros(nx+1,ny+2);
us = zeros(nx+1,ny+2); % "u-star"

% v (y-component velocity)
v = zeros(nx+2,ny+1);
vs = zeros(nx+2,ny+1); % "v-star"

% Pressure 
p = zeros(nx,ny);

% Values corresponding to pressure
pu = zeros(nx-1,ny);
pv = zeros(nx,ny-1);

% Storage of pressure field partial derivatives
dpx = zeros(nx-1,ny);
dpy = zeros(nx,ny-1);


Ru = zeros(nx,ny);
Rv = zeros(nx,ny);

% Centered u and v velocity vector components
uce = zeros(nx,ny);
vce = zeros(nx,ny);

uco = zeros(nx+1,ny+1);
vco = zeros(nx+1,ny+1);

% Passes geometric data to separate function which returns the toeplitz pressure
% matrix
if H == 1
    S = nx*ny; % Size of Laplacian
    L = Neumann_pressure_matrix(S,nx,ny,dxi2,dyi2);
    %L = pressure_matrix(nx,ny,dxi2,dyi2);
    L(1,:) = 0;
    L(1,1) = 1;

    % Invert the Laplacian (Pressure Matrix)
    LI = inv(L);
end


% [POSSIBLE INTITIAL CONDITION IMPLEMENTATION]

%% Plotting

% Mesh to plot velocity vectors on spatial domain
[Xce,Yce] = meshgrid(xce,yce);

% "Interp" - For vector visibility
%uvis_interp = interp2(Xce,Yce,uce,xx,yy,'cubic');
%vvis_interp = interp2(Xce,Yce,vce,xx,yy,'cubic');

[~,h_abs] = contourf(Xce',Yce',sqrt(uce.^2 + vce.^2)) ;

%hold on
%h_quiver = quiver(xx',yy',uvis_interp,vvis_interp,0.75,'Color',[1,1,1],'LineWidth',1);
xlim([0 Lx]); 
ylim([0 Ly]);

harrow = annotation('textarrow',[0.3 0.7],[0.96 0.96],"LineWidth",2);
haxes = gca;

%hold off

%% Outer Loop - Steps through "frames in time"  
for ii = 1:nt
%% Updating & Enforcing Boundary Conditions 

    % Compute Boundary Conditions (B.C.)
    u(:,1) = 2*ubcb - u(:,2);         % x-velocity bottom B.C.
    u(:,end) = 2*ubct - u(:,end-1);   % x-velocity top B.C.
    u(1,:) = 0;                       % x-velocity left B.C.
    u(end,:) = 0;                     % x-velocity right B.C.
    v(:,1) = 0;                       % y-velocity bottom B.C.
    v(:,end) = 0;                     % y-velocity top B.C.
    v(1,:) = 2*vbcl - v(2,:);         % y-velocity left B.C.
    v(end,:) = 2*vbcr - v(end-1,:);   % y-velocity right B.C.
   
    % Compute Boundary Conditions (B.C.) "Star" Matrices
    us(:,1) = 2*ubcb - u(:,2);         % x-velocity bottom B.C.
    us(:,end) = 2*ubct - u(:,end-1);   % x-velocity top B.C.
    us(1,:) = 0;                        % x-velocity left B.C.
    us(end,:) = 0;                      % x-velocity right B.C.
    vs(:,1) = 0;                        % y-velocity bottom B.C.
    vs(:,end) = 0;                      % y-velocity top B.C.
    vs(1,:) = 2*vbcl - v(2,:);         % y-velocity left B.C.
    vs(end,:) = 2*vbcr - v(end-1,:);   % y-velocity right B.C.
    % setting us and vs boundary conditions using u and v, rather than
    % using us and vs yields a solution that is exactly the same.
    % The only difference is that using u and v makes the ghost nodes in
    % this code match those in th Chinese code, so comparison becomes more
    % obvious
    % Computing the u and v at the center of the cell will yield more accurate results 
    uce = (u(1:end-1,2:end-1) + u(2:end,2:end-1))/2;
    vce = (v(2:end-1,1:end-1) + v(2:end-1,2:end))/2;
    
    for j = 1:(ny+1)
        for i = 1:(nx+1)
            uco(i,j) = (u(i,j)+u(i,j+1))*0.5; % u-values at the corners
        end
    end
    
    for i = 1:(nx+1)
        for j = 1:(ny+1)
            vco(i,j) = (v(i,j)+v(i+1,j))*0.5; % v values at the corners
        end
    end
    
%% "Predictor Step" - Spatial Data
 
    % Solve for u-star in the absence of pressure
    % This evalutation will not satisfy continuity 
    for j = 2:(ny+1)
        for i = 2:nx
            
            v_here = 0.25*(v(i+1,j) + v(i+1,j-1) + v(i,j-1) + v(i,j));
            
            % [FULL MOMENTUM EQ FOR REFERENCE] (P.6 OWKES)
            %us(i,j) = u(i,j) + dt*(nu*((u(i-1,j)-2*u(i,j)+u(i+1,j))*dxi2 + (u(i,j-1)-2*u(i,j)+u(i,j+1))*dyi2) ... 
            %- (u(i,j)*(u(i+1,j)-u(i-1,j))*0.5*dxi + v_here*(u(i,j+1)-u(i,j-1))*0.5*dyi)); 
            
            %% Differencing of Viscous and Convective Terms
            
            % Viscous Terms
            Lux = (u(i-1,j)-2*u(i,j)+u(i+1,j))*dxi2;
            Luy = (u(i,j-1)-2*u(i,j)+u(i,j+1))*dyi2;
            
            % Convective Terms
            Nux = u(i,j)*(uce(i,j-1)-uce(i-1,j-1))*dxi; % need the j-1 to keep indexing properly working
            Nuy = v_here*(uco(i-1,j)-uco(i-1,j-1))*dyi;
            
            % Momentum Conservation Finite Difference u-direction (EQ. 7 OWKES)
            us(i,j) = u(i,j) + dt*(nu*(Lux + Luy) - (Nux + Nuy)) ;
            
        end
    end
    
    
    
    % Solve for v-star in the absence of pressure
    % This evaluation will not satisfy continuity
    for j = 2:ny
        for i = 2:(nx+1)
            
            u_here = 0.25*(u(i,j+1) + u(i-1,j) + u(i,j) + u(i-1,j+1));
            
            % [FULL MOMENTUM EQ FOR REFERENCE] (P.6 OWKES)
            %vs(i,j) = v(i,j) + dt*(nu*((v(i-1,j)-2*v(i,j)+v(i+1,j))*dxi2 + (v(i,j-1)-2*v(i,j)+v(i,j+1))*dyi2) ... 
            %-(u_here*(v(i+1,j)-v(i-1,j))*0.5*dxi + v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*dyi)); 
             
            %% Differencing of Viscous and Convective Terms
            
            % Viscous Terms
            Lvx = (v(i-1,j)- 2*v(i,j) + v(i+1,j))*dxi2;
            Lvy = (v(i,j-1)- 2*v(i,j) + v(i,j+1))*dyi2;
            
            % Convective Terms
            Nvx = v(i,j)*(vce(i-1,j) - vce(i-1,j-1))*dxi; 
            Nvy = u_here*(vco(i,j-1) - vco(i-1,j-1))*dyi;
            
            % Momemntum Conservation Finite Difference v-direction (EQ. 12 OWKES)
            vs(i,j) = v(i,j) + dt*(nu*(Lvx + Lvy) - (Nvx + Nvy));
            
        end
    end
    
    
%% Computing div(u-star) - Center Difference (OPTION I)
    switch F
        case 1
            usce = (us(1:end-1,2:end-1)+us(2:end,2:end-1))/2;
            vsce = (vs(2:end-1,1:end-1)+vs(2:end-1,2:end))/2;

            % Partial u-star, u-component, with respect to x
            for j = 1:ny
                for i = 1:nx
                    if (i ~= nx)
                        Ru(i,j) = (usce(i+1,j)-usce(i,j))*dxi;
                    else
                        Ru(i,j) = (usce(i,j)-usce(i-1,j))*dxi;
                    end
                end
            end
    
            % Partial of u-star, ucomponent, with respect to y
            for j = 1:ny
                for i = 1:nx
                    if (j ~= ny)
                        Rv(i,j) = (vsce(i,j+1)-vsce(i,j))*dyi;
                    else
                        Rv(i,j) = (vsce(i,j)-vsce(i,j-1))*dyi;
                    end
                end
            end
            % END OPTION I
            
% Computing div(u-star) - Forward Difference (OPTION II) 
        case 2
            % u-component partial derviative (EQ 19 OWKES)
            for j = 2:(ny+1)
                for i = 1:nx
                    Ru(i,j-1) = (us(i+1,j) - us(i,j))*dxi; 
                end
            end
    
            % v-component partial derivative (EQ 19 OWKES)
            for j = 1:ny
                for i = 2:(nx+1)
                    Rv(i-1,j) = (vs(i,j+1) - vs(i,j))*dyi; 
                end
            end
            % END OPTION II
        otherwise
            disp('Return to line 23 and ensure that F is set to either 1 or 2')
    end

%% Corrector Step I - Pressure Contribution
   %  b = ((us(2:end,2:end-1)-us(1:end-1,2:end-1))/dx ...
       %+ (vs(2:end-1,2:end)-vs(2:end-1,1:end-1))/dy); 
    
    
    % Solve for p (using cosine transform, faster)
    %p = solvePoissonEquation_2dDCT(b,nx,ny,dx,dy);
    % Right Hand Side of Poisson Equation (EQ 6 OWKES)
    R = rho*dti*(Ru + Rv); % Summation of partial derviatives of component momentum functions
    if H == 1
        R = R(:);
        p = LI*R; % Inversion of Laplacian (EQ 17 OWKES)
        p = reshape(p,nx,ny); % Reformulation as matrix
    end
    
    % R takes care of rho and dt so I dont need them in the pressure
    % equation
    %u(2:end-1,2:end-1) = us(2:end-1,2:end-1) -  (p(2:end,:)-p(1:end-1,:))/dx;
    %v(2:end-1,2:end-1) = vs(2:end-1,2:end-1) -  (p(:,2:end)-p(:,1:end-1))/dy;
    % I think I need to make a separate loop for the edges of the pressure
    % expression
    
    % Conquer this task tomorrow
    % The 4 lines below properly populate the pressure values at the
    % corners of the domain
    
    % The method below kinda ends up at the same place, but the solution
    % process seems somewhat oscillatory
    
    % I think I need to fix the velocity vector plot in order to better
    % diagnose this problem, but I think I might have to iterate over the
    % pressure solutions or something?
    if H == 2
        for m = 1:150 % The iteration thing is DEFINITELY working..but why?
            p(1,1) = 1/(dx^2+dy^2)*(p(2,1)*dy^2 + p(1,2)*dx^2 - dx^2*dy^2*R(1,1));
            p(nx,1) = 1/(dx^2+dy^2)*(p(nx-1,1)*dy^2 + p(nx,2)*dx^2 - dx^2*dy^2*R(nx,1));
            p(1,ny) = 1/(dx^2+dy^2)*(p(2,ny)*dy^2 + p(1,ny-1)*dx^2 - dx^2*dy^2*R(1,ny));
            p(nx,ny) = 1/(dx^2+dy^2)*(p(nx-1,ny)*dy^2 + p(nx,ny-1)*dx^2 - dx^2*dy^2*R(nx,ny));
            for j = 2:ny-1
                p(1,j) = 1/(2*dx^2+dy^2)*(p(2,j)*dy^2 + (p(1,j+1)+p(1,j-1))*dx^2 - dx^2*dy^2*R(1,j)); % LHS
            end

            for j = 2:ny-1
                p(nx,j) = 1/(2*dx^2+dy^2)*(p(nx-1,j)*dy^2 + (p(nx,j+1)+p(nx,j-1))*dx^2 - dx^2*dy^2*R(nx,j)); % RHS
            end

            for i = 2:nx-1
                p(i,ny) = 1/(dx^2+2*dy^2)*(p(i,ny-1)*dx^2 + (p(i+1,ny)+p(i-1,ny))*dy^2 - dx^2*dy^2*R(i,ny)); %Top
            end

            for i = 2:nx-1
                p(i,1) = 1/(dx^2+2*dy^2)*(p(i,2)*dx^2 + (p(i+1,1)+p(i-1,1))*dy^2 - dx^2*dy^2*R(i,1)); %Bot
            end

            for i = 2:nx-1
                for j = 2:ny-1
                    p(i,j) = 1/(dx^2+dy^2)*((p(i,j+1)+p(i,j-1))*0.5*dx^2 + (p(i+1,j)+p(i-1,j))*0.5*dy^2 - dx^2*dy^2*R(i,j));
                end
            end
        end
    end
    % Ok, it kind of works...Looks bad, but definitely "works"
    % I will need to investigate more:
    % First, I will need to compare the matrix method with these settings
    % Doing so will tell me if this method is identical
    % Next, I need to see if iteration will be necessary
    
    
    % Solution of Pressure Poisson Equation
    %p = LI*R; % Inversion of Laplacian (EQ 17 OWKES)
    %p = reshape(p,nx,ny); % Reformulation as matrix
    
    % Pressure Location
    % Places the node in the center of the cell in the u and v direction
   
    for j = 1:ny
        for i = 1:(nx-1)
            pu(i,j) = (p(i,j) + p(i+1,j))*0.5;
        end
    end
    
    for j = 1:(ny-1)
        for i = 1:nx
            pv(i,j) = (p(i,j) + p(i,j+1))*0.5;
        end
    end
   
    % Finite Difference Computation of Pressure Gradient (EQ 21 & EQ 22
    % OWKES)
    % Computes pressure gradient for x and y components 
    % Accounts for vicinity to boundary
    for j = 1:ny
        for i = 1:(nx-1)
            if (i ~= (nx-1))
                dpx(i,j) = (pu(i+1,j)-pu(i,j))*dxi;
            else
                dpx(i,j) = (pu(i,j)-pu(i-1,j))*dxi;
            end
        end
    end
    
    for j = 1:(ny-1)
        for i = 1:nx
            if (j ~= (ny-1))
                dpy(i,j) = (pv(i,j+1)-pv(i,j))*dyi;
            else
                dpy(i,j) = (pv(i,j)-pv(i,j-1))*dyi;
            end
        end
    end
    
    %Accounting for pressure, the evaluation satisfies continuity. The
    %intermediate steps are used to compute u and v in Corrector Step II.

%% Corrector Step II - Computing u and v using u-star and v-star (P.9, OWKES)
% reintroduce the dt on top of rho when you finish comparing
% I originally removed it because the DCT pressure solution does not use it
   u(2:end-1,2:end-1) = us(2:end-1,2:end-1) - dt/rho*dpx(1:end,1:end); 
   v(2:end-1,2:end-1) = vs(2:end-1,2:end-1) - dt/rho*dpy(1:end,1:end);  

%% Cell-centering & Plotting Frame-by-Frame 
    for j = 1:ny
        for i = 1:nx
            uce(i,j) = 0.5*(u(i,j+1) + u(i+1,j+1));
            vce(i,j) = 0.5*(v(i+1,j+1) + v(i+1,j));
        end
    end
                  
%uvis = uce;
%vvis = vce;
  
%uvis_interp = interp2(Xce,Yce,uvis,xx,yy,'cubic');
%vvis_interp = interp2(Xce,Yce,vvis,xx,yy,'cubic');
  
h_abs.ZData = sqrt(uce.^2+vce.^2);

%h_quiver.UData = uvis_interp;
%h_quiver.VData = vvis_interp;
drawnow % Displays current contours

end