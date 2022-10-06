function L = pressure_matrix(nx,ny,dxi2,dyi2)

% nx and ny are the sizes of the block segments
% they are not the size of the resultant laplacian
% the size of the laplacian is nx^2
% so we just need to pass in nx and square root it
% obviously, this means that nx and ny need to
% be perfect squares
D = -2*(dxi2 + dyi2);
Dx = -(dxi2 + 2*dyi2);
Dy = -(2*dxi2 + dyi2);
Dxy = -(dxi2 + dyi2);
T1 = toeplitz([D dxi2 zeros(1,sqrt(nx)-2)]);
T2 = toeplitz([0 dyi2 zeros(1,sqrt(nx)-2)]);
I = eye(sqrt(nx),sqrt(ny));
L = kron(T1,I) + kron(I,T2);
L(1,1) = Dxy;
L(sqrt(nx),sqrt(ny)) = Dxy;
    for j = 2:(sqrt(nx)-1)
        for i = 2:(sqrt(ny)-1)
            if i==j
                L(i,j) = Dx;
            end
        end
    end
    L(end,end) = Dxy;
    L(nx-sqrt(nx)+1,ny-sqrt(ny)+1) = Dxy;
    for j = (nx-sqrt(nx)+2):(nx-1)
        for i = (ny-sqrt(ny)+2):(ny-1)
            if i==j
                L(i,j) = Dx;
            end
        end
    end

    for j = (nx-sqrt(nx)+2):(nx-1)
        for i = (ny-sqrt(ny)+2):(ny-1)
            if i==j
               L(i,j) = Dx;
            end
        end
    end

    c = 1;
    while (c <= (sqrt(nx)-2))
        L(c*sqrt(nx)+1,c*sqrt(ny)+1) = Dy;
        L(c*sqrt(nx)+sqrt(nx),c*sqrt(ny)+sqrt(ny)) = Dy;
        c = c+1;
    end
end

