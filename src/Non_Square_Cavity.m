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

% Calculation mode for Ru and Rv (Default is F = 1)
F = 1;

% Initial Boundary Conditions, B.C. (u and v direction)
ubcb = 0;
ubct = 1;
vbcl = 0;
vbcr = 0;

%% "Container" Space Discretization (xy-plane cross-section) - Indexing

% Number of cells
n = 100; 

% Number of cells in x and y axis
nx = n;
ny = n;

% Dimensions of spatial domain
Lx = 2;
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
nt = 300; 

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


Ru = zeros(nx,ny);
Rv = zeros(nx,ny);

% Centered u and v velocity vector components
uce = zeros(nx,ny);
vce = zeros(nx,ny);

uco = zeros(nx+1,ny+1);
vco = zeros(nx+1,ny+1);

% Passes geometric data to separate function which returns the toeplitz pressure
% matrix
L = pressure_matrix(nx,ny,dxi2,dyi2);

L(1,:) = 0;
L(1,1) = 1;
S = size(L);

% Invert the Laplacian (Pressure Matrix)
if (S(1) == S(2))
    LI = inv(L);
else
    LI = pinv(L);
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

hold off

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
    us(:,1) = 2*ubcb - us(:,2);         % x-velocity bottom B.C.
    us(:,end) = 2*ubct - us(:,end-1);   % x-velocity top B.C.
    us(1,:) = 0;                        % x-velocity left B.C.
    us(end,:) = 0;                      % x-velocity right B.C.
    vs(:,1) = 0;                        % y-velocity bottom B.C.
    vs(:,end) = 0;                      % y-velocity top B.C.
    vs(1,:) = 2*vbcl - vs(2,:);         % y-velocity left B.C.
    vs(end,:) = 2*vbcr - vs(end-1,:);   % y-velocity right B.C.
    
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
    
    % Right Hand Side of Poisson Equation (EQ 6 OWKES)
    R = rho*dti*(Ru + Rv); % Summation of partial derviatives of component momentum functions
    R = R(:);
    
    % Solution of Pressure Poisson Equation
    p = LI*R; % Inversion of Laplacian (EQ 17 OWKES)
    p = reshape(p,nx,ny); % Reformulation as matrix
    
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

   u(2:end-1,2:end-1) = us(2:end-1,2:end-1) - dt/rho*dpx(1:end,1:end); 
   v(2:end-1,2:end-1) = vs(2:end-1,2:end-1) - dt/rho*dpy(1:end,1:end);  

%% Cell-centering & Plotting Frame-by-Frame 
    for j = 1:ny
        for i = 1:nx
            uce(i,j) = 0.5*(u(i,j+1) + u(i+1,j+1));
            vce(i,j) = 0.5*(v(i+1,j+1) + v(i+1,j));
        end
    end
                  
uvis = uce;
vvis = vce;
  
%uvis_interp = interp2(Xce,Yce,uvis,xx,yy,'cubic');
%vvis_interp = interp2(Xce,Yce,vvis,xx,yy,'cubic');
  
h_abs.ZData = sqrt(uce.^2+vce.^2);

%h_quiver.UData = uvis_interp;
%h_quiver.VData = vvis_interp;
drawnow % Displays current contours

end