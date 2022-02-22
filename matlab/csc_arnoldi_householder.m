function [WH,b] = csc_arnoldi_householder(Av,Ar,Ac,v,m,LUv,pc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes de Arnoldi's Method with Householder
%
% Entries:
%   - Av,Ar,Ac : matrix A of dimensions n*n in CSC storage
%   - v : vector for compute Krylov subspace = span{v,Av,(A^2)v,...,(A^m-1)v}
%   - m : dimension of Krylov subspace
%   - LUv  : values of LU decomposition matrix (for SSOR or ILU(0)) in
%              CSC_packed storage, i.e., if pc==3 or pc==4
%            Or values of D or abs(D) if pc==1 or pc==2
%   - pc  : preconditioning type
%       - 0 : No preconditioning
%       - 1 : Diagonal preconditioning
%       - 2 : Absolute diagonal preconditioning
%       - 3 : Symmetric SOR (Symmetric Gauss-Seidel -> SSOR, w=1) or ILU(0)
%
% Output
%   - WH : Matrix with Hessemberg matrix in its upper+1 part, and W vectors
%   in its down part
%   - b : beta coefficients for computing Householder Projectors P = I - beta*w*w'
% 
% WH:  (n+1)x(m+1)
%    | h_10   h_11   h_12    .    .    .       h_1m  |
%    | 1.0    h_21   h_22    .    .    .       h_2m  |
%    | w_21   1.0    h_32    .    .    .       h_3m  |
%    | w_31   w_32   1.0     .    .    .       h_4m  |
%    |  .      .      .      .                  .    |
%    |  .      .      .           .             .    |
%    |  .      .      .                h_mm-1  h_mm  |
%    | w_m-11 w_m-12 w_m-13  .    .    1.0     h_m+1m|
%    | w_m1   w_m2   w_m3    .    .    w_mm+1  1.0   |
%    |  .      .      .      .    .     .      .     |
%    |  .      .      .      .    .     .      .     |
%    | w_n1   w_n2   w_n3    .    .     .     w_nm   |
%
% Note that w_ii=1.0, and shall not be stored, but that get computations
% quite difficult, instead the m w_ii are stored, which increase by m the
% real optimal storage.
%
% WH(1,1) = beta for GMRES residual computation
%
% Sergio Castiblanco
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NO PRECONDITIONING
if pc==0

n = length(v);
if m>n
    return
end
z = v;
WH = zeros(n+1,m+1);
b = zeros(1,m+1);
[WH(2:n+1,1),b(1)] = householderv(z);

if m<n
%Case when m<n H is computed completely out of if j<=m inside loop (line 57)

for j=1:m+1

    %Computing H
    WH(1:j-1,j) = z(1:j-1);
    WH(j,j) = z(j) - b(j)*(z(j:n)'*WH(j+1:n+1,j));
    
    if j<=m  %this if presumibly can be removed, if Q dimensions are (n,m+1)
        
        %Computing v_j, based in WH previous results
        % vj = P1*P2*...*PJ*ej
        % i.e. succesive outer products by right side
        v = zeros(n,1);
        q = v;
        q(j) = 1;
        for k = j:-1:1
            bqTw = b(k)*(q(k:n)'*WH(k+1:n+1,k));
            for i = k:n
                v(i) = q(i) - WH(i+1,k)*bqTw;
            end
            q = v;
        end

        %Computing z_j+1
        % z_j+1 = Pj*Pj-1*...*P2*P1*(A*v_j)
        % i.e succesive outer products by left side
        zaux = csc_matvec(Av,Ar,Ac,v);
        for k = 1:j
            bzTw = b(k)*(zaux(k:n)'*WH(k+1:n+1,k));
            for i = k:n
                z(i) = zaux(i) - WH(i+1,k)*bzTw;
            end
            zaux = z;
        end

        %Computing new Householder vector (W(j+1)) and new beta
        [WH(j+2:n+1,j+1),b(j+1)] = householderv(z(j+1:end));
    end
end

elseif m==n
%Case when m=n H is computed partially out of if j<=m inside loop (line 57)

for j=1:m+1

    %Computing H
    WH(1:j-1,j) = z(1:j-1);
        
    if j<=m  %this if presumibly can be removed, if Q dimensions are (n,m+1)
        
        WH(j,j) = z(j) - b(j)*(z(j:n)'*WH(j+1:n+1,j));
        %Computing v_j, stored in j_th column of Q
        % vj = P1*P2*...*PJ*ej
        % i.e. succesive outer products by right side
        v = zeros(n,1);
        q = v;
        q(j) = 1;
        for k = j:-1:1
            bqTw = b(k)*(q(k:n)'*WH(k+1:n+1,k));
            for i = k:n
                v(i) = q(i) - WH(i+1,k)*bqTw;
            end
            q = v;
        end

        %Computing z_j+1
        % z_j+1 = Pj*Pj-1*...*P2*P1*(A*v_j)
        % i.e succesive outer products by left side
        zaux = csc_matvec(Av,Ar,Ac,v);
        for k = 1:j
            bzTw = b(k)*(zaux(k:n)'*WH(k+1:n+1,k));
            for i = k:n
                z(i) = zaux(i) - WH(i+1,k)*bzTw;
            end
            zaux = z;
        end

        %Computing new Householder vector (W(j+1)) and new beta
        [WH(j+2:n+1,j+1),b(j+1)] = householderv(z(j+1:end));
    end
end

end

%DIAGONAL PRECONDITIONING or ABSOLUTE DIAGONAL PRECONDITIONING
% IN THIS CASE LUv = D, or LUv = abs(D) respectively
elseif pc==1 || pc==2

n = length(v);
if m>n
    return
end
z = v;
WH = zeros(n+1,m+1);
b = zeros(1,m+1);
[WH(2:n+1,1),b(1)] = householderv(z);

if m<n
%Case when m<n H is computed completely out of if j<=m inside loop (line 57)

for j=1:m+1

    %Computing H
    WH(1:j-1,j) = z(1:j-1);
    WH(j,j) = z(j) - b(j)*(z(j:n)'*WH(j+1:n+1,j));
    
    if j<=m  %this if presumibly can be removed, if Q dimensions are (n,m+1)
        
        %Computing v_j, based in WH previous results
        % vj = P1*P2*...*PJ*ej
        % i.e. succesive outer products by right side
        v = zeros(n,1);
        q = v;
        q(j) = 1;
        for k = j:-1:1
            bqTw = b(k)*(q(k:n)'*WH(k+1:n+1,k));
            for i = k:n
                v(i) = q(i) - WH(i+1,k)*bqTw;
            end
            q = v;
        end

        %Computing z_j+1
        % z_j+1 = Pj*Pj-1*...*P2*P1*(A*v_j)
        % i.e succesive outer products by left side
        zaux = v./LUv;
        zaux = csc_matvec(Av,Ar,Ac,zaux);
        for k = 1:j
            bzTw = b(k)*(zaux(k:n)'*WH(k+1:n+1,k));
            for i = k:n
                z(i) = zaux(i) - WH(i+1,k)*bzTw;
            end
            zaux = z;
        end

        %Computing new Householder vector (W(j+1)) and new beta
        [WH(j+2:n+1,j+1),b(j+1)] = householderv(z(j+1:end));
    end
end

elseif m==n
%Case when m=n H is computed partially out of if j<=m inside loop (line 57)

for j=1:m+1

    %Computing H
    WH(1:j-1,j) = z(1:j-1);
        
    if j<=m  %this if presumibly can be removed, if Q dimensions are (n,m+1)
        
        WH(j,j) = z(j) - b(j)*(z(j:n)'*WH(j+1:n+1,j));
        %Computing v_j, stored in j_th column of Q
        % vj = P1*P2*...*PJ*ej
        % i.e. succesive outer products by right side
        v = zeros(n,1);
        q = v;
        q(j) = 1;
        for k = j:-1:1
            bqTw = b(k)*(q(k:n)'*WH(k+1:n+1,k));
            for i = k:n
                v(i) = q(i) - WH(i+1,k)*bqTw;
            end
            q = v;
        end

        %Computing z_j+1
        % z_j+1 = Pj*Pj-1*...*P2*P1*(A*v_j)
        % i.e succesive outer products by left side
        zaux = v./LUv;
        zaux = csc_matvec(Av,Ar,Ac,zaux);
        for k = 1:j
            bzTw = b(k)*(zaux(k:n)'*WH(k+1:n+1,k));
            for i = k:n
                z(i) = zaux(i) - WH(i+1,k)*bzTw;
            end
            zaux = z;
        end

        %Computing new Householder vector (W(j+1)) and new beta
        [WH(j+2:n+1,j+1),b(j+1)] = householderv(z(j+1:end));
    end
end

end

%SSOR OR ILU(0) PRECONDITIONING
elseif pc==3 || pc==4

n = length(v);
if m>n
    return
end
z = v;
WH = zeros(n+1,m+1);
b = zeros(1,m+1);
[WH(2:n+1,1),b(1)] = householderv(z);

if m<n
%Case when m<n H is computed completely out of if j<=m inside loop (line 57)

for j=1:m+1

    %Computing H
    WH(1:j-1,j) = z(1:j-1);
    WH(j,j) = z(j) - b(j)*(z(j:n)'*WH(j+1:n+1,j));
    
    if j<=m  %this if presumibly can be removed, if Q dimensions are (n,m+1)
        
        %Computing v_j, based in WH previous results
        % vj = P1*P2*...*PJ*ej
        % i.e. succesive outer products by right side
        v = zeros(n,1);
        q = v;
        q(j) = 1;
        for k = j:-1:1
            bqTw = b(k)*(q(k:n)'*WH(k+1:n+1,k));
            for i = k:n
                v(i) = q(i) - WH(i+1,k)*bqTw;
            end
            q = v;
        end

        %Computing z_j+1
        % z_j+1 = Pj*Pj-1*...*P2*P1*(A*v_j)
        % i.e succesive outer products by left side
        zaux = csc_solpacklu(LUv, Ar, Ac, v);
        zaux = csc_matvec(Av,Ar,Ac,zaux);
        for k = 1:j
            bzTw = b(k)*(zaux(k:n)'*WH(k+1:n+1,k));
            for i = k:n
                z(i) = zaux(i) - WH(i+1,k)*bzTw;
            end
            zaux = z;
        end

        %Computing new Householder vector (W(j+1)) and new beta
        [WH(j+2:n+1,j+1),b(j+1)] = householderv(z(j+1:end));
    end
end

elseif m==n
%Case when m=n H is computed partially out of if j<=m inside loop (line 57)

for j=1:m+1

    %Computing H
    WH(1:j-1,j) = z(1:j-1);
        
    if j<=m  %this if presumibly can be removed, if Q dimensions are (n,m+1)
        
        WH(j,j) = z(j) - b(j)*(z(j:n)'*WH(j+1:n+1,j));
        %Computing v_j, stored in j_th column of Q
        % vj = P1*P2*...*PJ*ej
        % i.e. succesive outer products by right side
        v = zeros(n,1);
        q = v;
        q(j) = 1;
        for k = j:-1:1
            bqTw = b(k)*(q(k:n)'*WH(k+1:n+1,k));
            for i = k:n
                v(i) = q(i) - WH(i+1,k)*bqTw;
            end
            q = v;
        end

        %Computing z_j+1
        % z_j+1 = Pj*Pj-1*...*P2*P1*(A*v_j)
        % i.e succesive outer products by left side
        zaux = csc_solpacklu(LUv, Ar, Ac, v);
        zaux = csc_matvec(Av,Ar,Ac,zaux);
        for k = 1:j
            bzTw = b(k)*(zaux(k:n)'*WH(k+1:n+1,k));
            for i = k:n
                z(i) = zaux(i) - WH(i+1,k)*bzTw;
            end
            zaux = z;
        end

        %Computing new Householder vector (W(j+1)) and new beta
        [WH(j+2:n+1,j+1),b(j+1)] = householderv(z(j+1:end));
    end
end

end

end

end