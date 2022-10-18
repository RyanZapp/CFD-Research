clear
close all
clc
rho = 1 ;
nu = 0.008926 ;
nx = 50;
ny = 50;
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

nt = 50; % max number of timesteps

dt = dx*dy/(nu*(dx + dy));% timestep size (determined by CFL condition)
dt = 0.125*dxi2/nu; % This is a backup timestep size in case the other one
%fails
dt = 0.01;

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
% the things I did during this coding session were:
% I took the inversion of L out of the loop, so that it only had to happen
% once. I tried to use spectral decomposition to get L to invert without
% being singular... I think I accomplished that task, but now everything
% seems kinda worsse??? So chris and I will have to compare our laplacian
% to the laplacian developed in the chinese code to really see what is
% going on
%LI = inv(L); 
L(1,:) = 0;
L(1,1) = 1; % I think the reason that some codes do this is to combat the problematic eigenvalue at the end
[V,~] = eig(L);
%LI = inv(L);
lambda = eig(L);
lambda(end) = [];
lambda_reciprocal = 1./lambda;
lambda_reciprocal(end+1) = 0; % use a pseudo inverse and set (1/lambda)=0 for lambda=0 
D1_inv = diag(lambda_reciprocal);
LI = transpose(V)*D1_inv*V;
% basically, we have 1 eigenvector that is essentially 0, and that is our
% problem....Im not sure why it is zero, or if having a zero eigenvector is
% bad....all i know is that it is our current issue
 % we have p V D V^T
%p = V^T D^-1 V


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
     % R1 messes up on the third iteration. It is already big on the third iteration, which means us is big on the third iteration, which means
    % u must have been made big at the end of the 2nd iteration
    for j = 1:ny
        for i = 2:(nx+1)
            R2(i-1,j) = (vs(i,j+1) - vs(i,j))*0.5*dyi;
        end
    end
    % again, on the 3rd iteration, R3 is huge right from the get go, which means vs was huge when calculating it, which means
    % v must have been huge at the end of the 2nd iteration
    R = -rho/dt*(R1 + R2); % R actually becomes huge at the 2nd passthrough
    % R is so huge, so it does not matter if L is small, since R is going
    % to accelerate the values anyways
    R = R(:);
    %L(1,:) = 0;
    %L(1,1) = 1; % for some reason, this stopped the matrix from being singluar....now, the only issue is that, for some reason
    % the velocity for u is degenerating at the top boundary to zero, so
    % the cavity flow quickly comes to a halt
    % I do not know if the entries of R will match up with the
    % corresponding entries of the laplacian
    p = L\R;
    %p = LI*R; % Rather than using backslash, we can just invert the laplacian a single time before we ever enter the loop
    % when i do L\R, I get a better looking result than when I spectrally
    % decompose L and invert it that way
    % In addition, when I use inv(L), it messes everything up, presumably
    % because the "inv" command cannot handle near singular matrices...What
    % ends up happening is I get an all constant matrix, so "inv"
    % definitely is not the move...To double check, I should ask chris what
    % happens when he uses the inv command
    p = reshape(p,nx,ny); % p is huge at the 2nd go around, so thatt means 
   % LL = inv(L); % L stays constant. R is the only thing that does not stay constant, so we can probably just "pre-invert"L
    
    % after I fix L, I just need to add the pressure derivatives back into
    % the solution
    
    for j = 2:(ny+1)
        for i = 2:nx
            u(i,j) = us(i,j) - dt/rho*(p(i,j-1) - p(i-1,j-1))*0.5*dxi; % I think this one is correct, but it would not hurt to double check
        end % so pressure somehow becomes huge, which is what is influencing u to get huge, which in turn, influences us to get huge
       u(:,1) = 2*ubcb - u(:,2);
    end
    
    u(:,end) = 2*ubct - u(:,end-1);
 
    u(1,:) = 0;
    u(end,:) = 0;
    for j = 2:ny
        for i = 2:(nx+1)
            v(i,j) = vs(i,j) - dt/rho*(p(i-1,j) - p(i-1,j-1))*0.5*dyi; % I think this is right, but I should double check it
        end
         v(1,:) = 2*vbcl - v(2,:);         % y-velocity left B.C.
         v(end,:) = 2*vbcr - v(end-1,:);   % y-velocity right B.C.
    end
    v(:,1) = 0;                     
    v(:,end) = 0;  
    % something weird is happening in the middle row of the v matrix
    for j = 1:ny
        for i = 1:nx
            uce(i,j) = 0.5*(u(i,j+1) + u(i+1,j+1));
            vce(i,j) = 0.5*(v(i+1,j+1) + v(i+1,j));
        end
       
    end
                  
  
    h_abs.ZData = sqrt(uce.^2+vce.^2);
    drawnow
    pause(0.025)
% There are more efficient ways to do this code, but let us just make sure
% it works for now

% The pressure matrix has issues
% I am not sure why for the v matrix, we initialized the Neumann
% boundary conditions in the loop...but in the u matrix, we did the bottom
% Neumann B.C. in the loop and the top Neumann B.C. out of the loop

% I honestly think that both can go outside of the loop



% Bottom line is, the velocity values are exploding probably due to a
% singularity sometwhere)

% The changes I made were on lines 197 and line 208. I realized that I was
% not multiplying the pressure term by dt. I was looking at videos and at
% the papaer, and apparently, that is what you are supposed to do...not
% 100% sure why because it did not 100% come up in the derivation.


% The reason that things were (and possibly still are) exploding to
% infinity is because when we assemble the R matrix, we divide by
% dt...because dt was SUPER small, it was making R super huge. although
% L^-1 was small, R was big enough to the point where it did not matter and
% L^-1*R was really big... and because p=L^-1*R, obviously, p is super big
% as a result. And then becasue we were adding the derivative of p to us,
% even though us was correct (i.e. it was small), obviously, adding a huge
% value in the form of p, to a small value in the form of us, still
% generates a huge value.

% That is when I looked at all of the formulas and things onliner and I
% realized that when we add the pressure derivative back into us and vs, we
% are supposed to multiply by dt to offset the absolute massive value of p.

% doing this actually ended up working, however, it sort of just staved off the
% inevitable, and things would run, but still end up crashing later

% Then I tried redoing all the combinations of the laplacian to see if that
% would help... WHat I noticed (and I will indeed double check this) is
% that using my own inverse, or pseudo inverse, or anything for that
% matter, would lead to a contour plot that wasnt even close to right. So I
% started going back to using the backslash command. That seemed to produce
% a musch more accurate plot...hoiwever, the matrix was still badly
% conditioned. IN order to combat that, I treid to use the other piece of
% information from the other paper where we set the top row of the
% laplacian equal to all zeros and then set the first entry in the matrix
% equal to 1 in order to stop the zero eigenvalue from appearing. THis got
% rid og the warning messages that the matrix was singular, which was nice,
% but the solution just appeared to oscillate through the first stages of
% the process..

% Now that I am looking back on it, I think what might be happening is the
% fact that maybe the pressure matrix is getting bigger and bigger until we
% reach the point where we just get NaNs, so no matter what, even if the
% computations are correct, once we exceed a certain value size on the
% computer, we will jsut get NaN computations or inf computations


% WHat this means is that we can probably put a break point just after line
% 171 to investigate how line 171 (which is the initialization of R)
% changes as the simulation continues....this might be able to tell us if
% we are running into an issue similar to what happened at the lab, where
% the numbers were just too big for t


% A final note on the matrix choice: I think that the backslash command or
% the linsolve commands are our best bets becasue my inversion method is
% trying to output complex values when the  grid gets
% refined....additionally, the plots for my inversion method look a bit
% bogus, but when I use the backslash command, they look like they are on
% the right track, so we should stick with that for now


% other than checking the increase in size of R, the only other thing that
% needs to be checked is how our outputs for u,v,uce, and vce compare to
% the outputs from the chinese github


% The only difference between when I set the values in the first row of the
% laplacian vs when I dont is that when I do set them, the matrix is not
% badly conditioned.......however, the plots between the two look almost
% identical
end