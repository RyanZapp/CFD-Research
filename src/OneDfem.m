clear
close all
clc


% Try using a ton of grid points, or try removing the dx
% We are solving u'' + u = f , u(a) = 0 , u(b) = 0
% Initialize parameters
a = 0;
b = 1;
N = 10000;
x = linspace(a,b,N)';
dx = (b-a)/(N-1);
S = length(x);
A = zeros(S,S);
% Initialize forcing function
f = x;
A(1,1) = 1; % Initialize the boundary shape function
A(end,end) = 1; % Initialize boundary shape function
for i = 2:S-1
	A(i-1,i) = 0; % This line is not necessary
	A(i,i) = 1;
	A(i+1) = 0; % This line is not necessary
end

% The lines labeled not necessary are just for the 
% sake of helping with conceptual understanding

% Compute the derivatives of the shape functions
B = zeros(S-1,S-1);
B = diff(A,1,1); % computes a forward difference down the columns of A
B(end+1,:) = B(end,:); % Apply a backward difference

% At this point, the columns of A represent the shape
% functions
M = zeros(S,1);
% Compute <f.phi>
for i = 1:S
	M(i) = sum(f.*A(:,i));
end

% Initialize Gram matrix
G = zeros(S,S);

for j = 1:S
	for i = 1:S
        if i == j
            G(j,i) = 2/3*dx + 2/dx;
        end
        if (i == j+1) || (i == j-1)
            G(j,i) = dx/6 - 1/dx;
        end
                
		%G(j,i) = sum((A(:,i).*A(:,j) - B(:,i).*B(:,j)));
	end
end
%G(1,:) = 0; % reenforce the Dirichlet B.C. on the test function;
%G(2,1) = 0;
%G(end,:) = 0; % reenforce the Dirichlet B.C. on the test function;
%G(end-1,end) = 0;
Ginv = pinv(G);
disp('Condition Number of G:')
disp(cond(G))
% [V,D] = eig(G);
% D = D(:);
% D(1333:1335) = 0.1;
% D = reshape(D,S,S);
% G = V*D*transpose(V);
%c = Ginv*M;
c = G\M;

u = zeros(S,1);
u_temp = zeros(S,1);
for i = 1:S
    u_temp = c(i)*A(:,i);
	u = u_temp + u;
end
%u_analytical = x - sin(x)/sin(1);
u_analytical = cos(x) - cot(1)*sin(x) + x;
figure
plot(x,u_analytical,'r')
hold on
plot(x,u,'bo')
xlabel('x')
ylabel('u(x)')
title('Numerical vs Analytical Solution')
legend('Analytical','Numerical','location','B')
grid on
% I might have to do a symbolic integration between f and phi as well