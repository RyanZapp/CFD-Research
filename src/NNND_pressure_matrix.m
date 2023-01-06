function L = NNND_pressure_matrix(nx,ny,dxi2,dyi2)

D = -2*(dxi2 + dyi2); % Is used in both
DxN = -(dxi2 + 2*dyi2); % N stands for Neumann
DyN = -(2*dxi2 + dyi2); % N stands for Neumann
DxyN = -(dxi2 + dyi2); % For the Neumann condition in the upper left block
DxyD = -(2*dxi2 + dyi2); % For the Dirichlet condition in the lower left block


T1 = toeplitz([D dxi2 zeros(1,nx-2)]);
T2 = toeplitz([0 dyi2 zeros(1,nx-2)]);
I = eye(nx,ny);
L = kron(T1,I) + kron(I,T2);

L(1,1) = DxyN;
L(nx,ny) = DxyN;
    for j = 2:(nx-1)
        for i = 2:(ny-1)
            if i==j
                L(i,j) = DxN;
            end
        end
    end
L(end,end) = DxyD;
L(nx^2-nx+1,ny^2-ny+1) = DxyD;
    for j = (nx^2-nx+2):(nx^2-1)
        for i = (ny^2-ny+2):(ny^2-1)
            if i==j
                L(i,j) = D;
            end
        end
    end
 

c = 1;
    while (c <= (nx-2))
        L(c*nx+1,c*ny+1) = DyN;
        L(c*nx+nx,c*ny+ny) = DyN;
        c = c+1;
    end
