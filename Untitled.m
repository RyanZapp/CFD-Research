clear
close all
clc
rho = 1 ;
nu = 0.008926 ;
n = 150; % while we work with square regions, lets just set this for the sake of quick testing
nx = n;
ny = n;
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

nt = 1500; % max number of timesteps

%dt = dx*dy/(nu*(dx + dy));% timestep size (determined by CFL condition)
%dt = 0.125*dxi2/nu; % This is a backup timestep size in case the other one
%fails
dt = 0.001;

% Specify location of cell centers (this is where we will be plotting)
xce = ((1:nx)-0.5)*dx; % length(xce) = nx
yce = ((1:ny)-0.5)*dy; % length(yce) = ny

% Create vector that houses spatial x coordinates (starting from bottom
% left corner of physical domain)
xco = (0:nx)*dx;  % length(xco) = nx + 1
yco = (0:ny)*dy; % length(yco) = ny + 1
delta = 0.05;
[xx,yy] = meshgrid(0:delta:Lx,0:delta:Ly);
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
uco = zeros(nx+1,ny+1);

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
% The 2 lines implemented below are implemented in a lot of the codes that
% I see
L(1,:) = 0;
L(1,1) = 1; % I think the reason that some codes do this is to combat the problematic eigenvalue at the end
% [V,~] = eig(L);
% %LI = inv(L);
% lambda = eig(L);
% lambda(end) = [];
% lambda_reciprocal = 1./lambda;
% lambda_reciprocal(end+1) = 0; % use a pseudo inverse and set (1/lambda)=0 for lambda=0 
% D1_inv = diag(lambda_reciprocal);
% LI = transpose(V)*D1_inv*V;
% basically, we have 1 eigenvector that is essentially 0, and that is our
% problem....Im not sure why it is zero, or if having a zero eigenvector is
% bad....all i know is that it is our current issue
 % we have p V D V^T
%p = V^T D^-1 V

%uuu = zeros(nx,ny);
% to access the interior points of each matrix, type the following:
% u(2:end-1,2:end-1);
% v(2:end-1,2:end-1);
% p(1:end,1:end);

% at this point, if you have a nonzero initial condition, you would
% populate u and v with that initial condition here

% however, our initial condition is the zeros matrix, so we will proceed.

% Create the meshgrid that you will eventually plot over
[Xce,Yce] = meshgrid(xce,yce);
uvis_interp = interp2(Xce,Yce,uce,xx,yy,'cubic');
vvis_interp = interp2(Xce,Yce,vce,xx,yy,'cubic');
%figure
[~,h_abs] = contourf(Xce',Yce',sqrt(uce.^2+vce.^2));
hold on
h_quiver = quiver(xx',yy',uvis_interp,vvis_interp,0.75,'Color',[1,1,1],'LineWidth',1);
xlim([0 Lx]); 
ylim([0 Ly]);
harrow = annotation('textarrow',[0.3 0.7],[0.96 0.96],"LineWidth",2);
haxes = gca;
%haxes.XTick = [];
%haxes.YTick = [];
hold off

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
    % No differences up to this point
    % These matrices are exteremely messed up********
    % note that us is somehow coming out wrong
    % Chris and I need to fix them before proceeding
    %uuu = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2;
    for j = 2:(ny-1)
        for i = 1:(nx-1)
        uce(i,j-1) = (u(i,j)+u(i+1,j))*0.5; % uce uses the exact u values, so its always accuracte
        end
    end
    vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2;
    for j = 1:(ny+1)
        for i = 1:(nx+1)
            uco(i,j) = (u(i,j)+u(i,j+1))*0.5; % u values at the corners
        end
    end
    % solve for u* (u* does not satisfy the continuity equation)
    for j = 2:(ny+1)
        for i = 2:nx
            v_here = 0.25*(v(i+1,j) + v(i+1,j-1) + v(i,j-1) + v(i,j));
            us(i,j) = u(i,j) + dt*(nu*((u(i-1,j)-2*u(i,j)+u(i+1,j))*dxi2 + (u(i,j-1)-2*u(i,j)+u(i,j+1))*dyi2) ... 
            - (u(i,j)*(u(i+1,j)-u(i-1,j))*0.5*dxi + v_here*(u(i,j+1)-u(i,j-1))*0.5*dyi));
            Lux(i-1,j-1) = (u(i-1,j)-2*u(i,j)+u(i+1,j))*dxi2;
            Luy(i-1,j-1) = (u(i,j-1)-2*u(i,j)+u(i,j+1))*dyi2;
            Nux(i-1,j-1) = u(i,j)*(uce(i,j-1)-uce(i-1,j-1))*dxi; % need the j-1 to keep indexing properly working
            Nuy(i-1,j-1) = v_here*(uco(i-1,j)-uco(i-1,j-1))*dyi;
            Lux1 = (u(i-1,j)-2*u(i,j)+u(i+1,j))*dxi2;
            Luy1 = (u(i,j-1)-2*u(i,j)+u(i,j+1))*dyi2;
            Nux1 = u(i,j)*(uce(i,j-1)-uce(i-1,j-1))*dxi; % need the j-1 to keep indexing properly working
            Nuy1 = v_here*(uco(i-1,j)-uco(i-1,j-1))*dyi;
            %us(i,j) = u(i,j) + dt*(nu*(Lux1+Luy1) - (Nux1+Nuy1)); % our us does differ from the chinese output of us
            % I checked it and its literally just because we are
            % multiplying by nu, which means that our formulation for
            % v_here is sufficient
            
            
            % lets keep Lux, Luy, Nux, and Nuy separate for now, but when
            % the time comes, we should combine them into a matrix called
            % du for the sake of computational efficiency
            
            
            % or maybe we should keep them alone, but as scalars, and then
            % just add them to a matrix called du
            
            
            
            % I might as well add all 4 of these terms and store them in a
            % single quantity
            % The only thing I can change about this loop without putting
            % things into convective form is to create vco and then use it
            % to create v_here
        end
        %us(:,1) = 2*ubcb - us(:,2);         % x-velocity bottom B.C.
    end
    %us(:,end) = 2*ubct - us(:,end-1);   % x-velocity top B.C.
    %us(1,:) = 0;                       % x-velocity left B.C.
    %us(end,:) = 0;                     % x-velocity right B.C.
    % we put the dirichlet B.C. again becasue the Neumann puts the velocity
    % at the top corners of the box as being equal to the lid velocity,
    % when in reality, we should have a no-slip at the corners.
    % Reintroducing the Dirichlet B.C. fixes this issue
    
%     Lu = nu*dt*(Lux+Luy);
%     Nu = dt*(Nux+Nuy);
%     u(2:end-1,2:end-1) = u(2:end-1,2:end-1) + Lu - Nu; % this matches my modified us computation
    % the only thing that "changes things" is when I remodify u later down
    % the line0
    
    % solve for v* (u* does not satisfy the continuity equation)
    for j = 2:ny
        for i = 2:(nx+1)
            u_here = 0.25*(u(i,j+1) + u(i-1,j) + u(i,j) + u(i-1,j+1));
            vs(i,j) = v(i,j) + dt*(nu*((v(i-1,j)-2*v(i,j)+v(i+1,j))*dxi2 + (v(i,j-1)-2*v(i,j)+v(i,j+1))*dyi2) ... 
            -(u_here*(v(i+1,j)-v(i-1,j))*0.5*dxi + v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*dyi));
            Lvx(i-1,j-1) = (v(i-1,j)-2*v(i,j)+v(i+1,j))*dxi2;
            Lvy(i-1,j-1) = (v(i,j-1)-2*v(i,j)+v(i,j+1))*dyi2;
        end
       % vs(1,:) = 2*vbcl - vs(2,:);         % y-velocity left B.C.
       % vs(end,:) = 2*vbcr - vs(end-1,:);   % y-velocity right B.C.
    end
   % vs(:,1) = 0;                       % y-velocity bottom B.C.
   % vs(:,end) = 0;                     % y-velocity top B.C.
    Lv = nu*dt*(Lvx+Lvy);
    usce = (us(1:end-1,2:end-1)+us(2:end,2:end-1))/2;
    vsce = (vs(2:end-1,1:end-1)+vs(2:end-1,2:end))/2;
%     for j = 2:(ny+1)
%         for i = 1:nx
%             R1(i,j-1) = (us(i+1,j) - us(i,j))*dxi; % I think we only need to scale by dxi because the
%             % points are only a distance of dx away from each other
%         end
%     end
    
    for j = 1:ny
        for i = 1:nx
            if (i ~= nx)
                R1(i,j) = (usce(i+1,j)-usce(i,j))*dxi;
            else
                R1(i,j) = (usce(i,j)-usce(i-1,j))*dxi;
            end
        end
    end
     % R1 messes up on the third iteration. It is already big on the third iteration, which means us is big on the third iteration, which means
    % u must have been made big at the end of the 2nd iteration
    
%     for j = 1:ny
%         for i = 2:(nx+1)
%             R2(i-1,j) = (vs(i,j+1) - vs(i,j))*dyi; % we want to get values at the cell centers
%         end
%     end
    
    for j = 1:ny
        for i = 1:nx
            if (j ~= ny)
                R2(i,j) = (vsce(i,j+1)-vsce(i,j))*dyi;
            else
                R2(i,j) = (vsce(i,j)-vsce(i,j-1))*dyi;
            end
        end
    end
    % again, on the 3rd iteration, R3 is huge right from the get go, which means vs was huge when calculating it, which means
    % v must have been huge at the end of the 2nd iteration
    R = -rho/dt*(R1 + R2); % R actually becomes huge at the 2nd passthrough
    % R is so huge, so it does not matter if L is small, since R is going
    % to accelerate the values anyways
    R = R(:);
    %R(1) = 0;
    %L(1,:) = 0;
    %L(1,1) = 1; % for some reason, this stopped the matrix from being singluar....now, the only issue is that, for some reason
    % the velocity for u is degenerating at the top boundary to zero, so
    % the cavity flow quickly comes to a halt
    % I do not know if the entries of R will match up with the
    % corresponding entries of the laplacian
    b = ((us(2:end,2:end-1)-us(1:end-1,2:end-1))/dx ...
        + (vs(2:end-1,2:end)-vs(2:end-1,1:end-1))/dy); % this is just taking the sum of partial derivatives as the divergence
    % why dont they scale b by rho and dt???
    
    % Solve for p (using cosine transform, faster)
    p = solvePoissonEquation_2dDCT(b,nx,ny,dx,dy); % this 
   % p = L\R;
    %p = LI*R; % Rather than using backslash, we can just invert the laplacian a single time before we ever enter the loop
    % when i do L\R, I get a better looking result than when I spectrally
    % decompose L and invert it that way
    % In addition, when I use inv(L), it messes everything up, presumably
    % because the "inv" command cannot handle near singular matrices...What
    % ends up happening is I get an all constant matrix, so "inv"
    % definitely is not the move...To double check, I should ask chris what
    % happens when he uses the inv command
    %p = reshape(p,nx,ny); % p is huge at the 2nd go around, so thatt means 
   % LL = inv(L); % L stays constant. R is the only thing that does not stay constant, so we can probably just "pre-invert"L
   % p(1,1) = 0;
    % after I fix L, I just need to add the pressure derivatives back into
    % the solution
%     for j = 1:ny
%         for i = 1:(nx-1)
%             pu(i,j) = (p(i,j) + p(i+1,j))*0.5;
%         end
%     end
%     for j = 1:(ny-1)
%         for i = 1:nx
%             pv(i,j) = (p(i,j) + p(i,j+1))*0.5;
%         end
%     end
    for j = 2:(ny+1)
        for i = 2:nx
            u(i,j) = us(i,j) - 1/rho*(p(i,j-1) - p(i-1,j-1))*dxi; % I think this one is correct, but it would not hurt to double check
        end % so pressure somehow becomes huge, which is what is influencing u to get huge, which in turn, influences us to get huge
       %u(:,1) = 2*ubcb - u(:,2);
    end
%     for j = 1:ny
%         for i = 1:(nx-1)
%             if (i ~= (nx-1))
%                 dpx(i,j) = (pu(i+1,j)-pu(i,j))*dxi;
%             else
%                 dpx(i,j) = (pu(i,j)-pu(i-1,j))*dxi;
%             end
%         end
%     end
%     
%     for j = 1:(ny-1)
%         for i = 1:nx
%             if (j ~= (ny-1))
%                 dpy(i,j) = (pv(i,j+1)-pv(i,j))*dyi;
%             else
%                 dpy(i,j) = (pv(i,j)-pv(i,j-1))*dyi;
%             end
%         end
%     end
    %u(:,end) = 2*ubct - u(:,end-1);
 
    %u(1,:) = 0;
    %u(end,:) = 0;
%   u(2:end-1,2:end-1) = us(2:end-1,2:end-1) - 1/rho*dpx(1:end,1:end);
%   v(2:end-1,2:end-1) = vs(2:end-1,2:end-1) - 1/rho*dpy(1:end,1:end);
%   
%     for j = 2:(ny+1)
%         for i = 2:nx
%             u(i,j) = us(i,j) - dt/rho*dpx(i-1,j);
%         end
%     end
%     
%     for j = 2:ny
%         for i = 2:(nx+1)
%             v(i,j) = vs(i,j) - dt/rho*dpy(i,j);
%         end
%     end
    
    for j = 2:ny
        for i = 2:(nx+1)
            % I think there should be a dt multiplying the pressure
            % derivative
            v(i,j) = vs(i,j) - 1/rho*(p(i-1,j) - p(i-1,j-1))*dyi; %I only need to divide by dy because the differencing 
            % points are only dy apart, not 2dy apart (I was previously
            % using 2dy because I was robotically using the central 
        end
         v(1,:) = 2*vbcl - v(2,:);         % y-velocity left B.C.
         v(end,:) = 2*vbcr - v(end-1,:);   % y-velocity right B.C.
    end
%     v(:,1) = 0;                     
%     v(:,end) = 0;  
    % something weird is happening in the middle row of the v matrix
    for j = 1:ny
        for i = 1:nx
            uce(i,j) = 0.5*(u(i,j+1) + u(i+1,j+1));
            vce(i,j) = 0.5*(v(i+1,j+1) + v(i+1,j));
        end
       
    end
                  
  uvis = uce;
  vvis = vce;
  %uvis(:,1:ceil(end/2)) = 0;
  %vvis(:,1:ceil(end/2)) = 0;
  uvis_interp = interp2(Xce,Yce,uvis,xx,yy,'cubic');
  vvis_interp = interp2(Xce,Yce,vvis,xx,yy,'cubic');
    h_abs.ZData = sqrt(uce.^2+vce.^2);
   % quiver(Xce',Yce',uce,vce,3,'Color',[1,1,1]);
    %hold on
    h_quiver.UData = uvis_interp;
    h_quiver.VData = vvis_interp;
    drawnow
    %pause(0.0000000000000001)
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




% Big interesting update....at least after the first iteration, it seems
% that the output of every version of us that I have is equivalent to the
% intermediate u velocity obtained in the chinese code (as long as I
% replace multiplication by nu with division by Re, where Re = 300)

% This kind of scares me and leads me to believe that the pressure step may
% be cracked?

% so maybe the best move is to investigate the output of the pressure
% derivatives and compare them to that of the chinese code to see if things
% match up initially??

% also, I should probably formulate a new vs just so things remain coherent
% and easy to test






% Literally, all you need to do to answer your pressing question is to swap
% out your pressure solver for the solve poisson equation and then if
% things end up looking great, then you know that the poisson solver is the
% issue


% then you can work on cleaning some of the code up using the strategies
% you learned, or you can keep it as is, since all you would have to do is
% compare your us and vs outputs at each step to that of the chinese code
% to figure out once and for all if all of these "corner" and "center
% velocity" obtaining schemes are equivalent


% I tried implementing the chinese pressure solver in my code and it went
% smoothly, but the output pressure matrix did not match up

% I checked us and it matches with the u in the chinese code that is being


% takeaways, increasing nu means that we need to have a smaller dt value
% (this means that dt depends on nu)

% additionally, I havent checked for 1-1 correspondance, between the
% chinese equations and our equations when we use the reynolds number
% formulation of the equations and we just use their pressure matrix,
% butthe solutions seem to match]



% the for sure thing is that using vs not using the non dimensional navier
% stokes equations does not matter. Additionally, using vs not using the
% convective form does not seem to matter in this scenario (as I suspected)

% we will still have to check, but also, using the corner values and center
% values vs our original method does not seem to matter either (this needs
% to be heavily tested and verified though)


% the only thing that appears to be wrong is that our poisson solver was
% booty....which does not make much sense because I am fiarly confident
% that I assembled the laplacian properly, and I am effectively inverting
% it)


% This makes me want to try and use the uce values to calculate us at the
% center of the cells, because finite differencing that way makes a lot
% more sense to me

% additionally, I am not sure why we had to remove dt from the
% reintroduction of pressure into the equation....is it a special thing due
% to the way that the chinese code solved the pressure, or is the multiple
% of dt actually not supposed to be there? (we will need to consult the
% literature to answer that question)

% additionally, when reintroducing pressure into the mix, I honestly think
% that we should obtain the pressure at the cell wallls and then compute
% the derivative so that we can m,athematically say that it is bering
% introduced at the right spt

% That being said, if the chinese method of computing us and our method of
% computing us are equivalent, then that means my central difference
% between the cell walls is the same as computing the value at the cell
% wall by differencing the values at the left and right cell centers of the
% cell wall, which kind of makes sense, but also kind of does not make
% sense



% more testing will be required to understand this



% I tried my ideas to fix the pressure solver and they dont work


% we should clean up the code and remove the BC where we dont need them

% then we should stash the pressure loop ideas somewhere else

% then, things should be cleaned up enough for us to begin trying to
% understand why our pressure solver does not work
b = ((u(2:end,2:end-1)-u(1:end-1,2:end-1))/dx ...
    + (v(2:end-1,2:end)-v(2:end-1,1:end-1))/dy);
if (b>=1e-3)
    disp("Warning, divergence is nonzero!")
    disp("The error occured during the following iteration:")
    disp(ii)
end

end