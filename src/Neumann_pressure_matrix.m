function L = Neumann_pressure_matrix(S,nx,ny,dxi2,dyi2)
    % Initialize parameters
    D = -2*(dxi2 + dyi2);
    Dx = -(dxi2 + 2*dyi2);
    Dy = -(2*dxi2 + dyi2);
    Dxy = -(dxi2 + dyi2);
    m = 1;
    L = zeros(S,S);
    % Assemble first block of matrix
    V1 = zeros(1,ny);
    V1(1) = Dx;
    V1(2) = dyi2;
    D1 = toeplitz(V1);
    D1(1,1) = Dxy;
    D1(end,end) = Dxy;
    % Assemble second block of matrix
    V2 = zeros(1,ny);
    V2(1) = D;
    V2(2) = dyi2;
    D2 = toeplitz(V2);
    D2(1,1) = Dy;
    D2(end,end) = Dy;
    % Assemble third block of matrix
    V3 = zeros(1,ny);
    V3(1) = dxi2;
    D3 = toeplitz(V3);
    % Assemble fourth block of matrix
    D4 = zeros(ny,ny);
    % Assemble The rows of the block matrix
    R1 = [D1 D3];
    Rend = [D3 D1];
    Rcent = [D3 D2 D3];
    S1 = size(R1);
    Send = size(Rend);
    Scent = size(Rcent);
    % Insert first and last block rows into Laplacian
    L(1:S1(1),1:S1(2)) = R1;
    L(end-Send(1)+1:end,end-Send(2)+1:end) = Rend;
    % Insert remaining block rows into Laplacian
    K = nx-3;
    if K ~= 0
        for i = 1:K
            Rcent = horzcat(Rcent,D4);
        end
    end
    L(S1(1)+1:S1(1)+Scent(1),1:end) = Rcent;
    while m <= K
        Rcent = circshift(Rcent,ny,2);
        L(S1(1)+m*Scent(1)+1:S1(1)+(m+1)*Scent(1),1:end) = Rcent;
        m = m + 1;
    end
end