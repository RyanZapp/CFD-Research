function L = pressure_matrix(nx,ny,dxi2,dyi2)

D = -2*(dxi2 + dyi2);
Dx = -(dxi2 + 2*dyi2);
Dy = -(2*dxi2 + dyi2);
Dxy = -(dxi2 + dyi2);
T1 = toeplitz([D dxi2 zeros(1,nx-2)]);
T2 = toeplitz([0 dyi2 zeros(1,nx-2)]);
I = eye(nx,ny);
L = kron(T1,I) + kron(I,T2);

L(1,1) = Dxy;
L(nx,ny) = Dxy;
    for j = 2:(nx-1)
        for i = 2:(ny-1)
            if i==j
                L(i,j) = Dx;
            end
        end
    end
L(end,end) = Dxy;
L(nx^2-nx+1,ny^2-ny+1) = Dxy;
    for j = (nx^2-nx+2):(nx^2-1)
        for i = (ny^2-ny+2):(ny^2-1)
            if i==j
                L(i,j) = Dx;
            end
        end
    end
 

c = 1;
    while (c <= (nx-2))
        L(c*nx+1,c*ny+1) = Dy;
        L(c*nx+nx,c*ny+ny) = Dy;
        c = c+1;
    end
