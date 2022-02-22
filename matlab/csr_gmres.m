function [x,t,res] = csr_gmres(Av,Ac,Ar,b,x,m,niter,tol,LUv,pc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function solves the system Ax=b with the Restarted Generalized
% Minimal Residual (Restarted GMRES)
%
% Entries:
%     Av,Ac,Ar : Matrix A in CSR storage
%     b : Right hand side vector
%     x : First guest for the solution
%     m : Krylov space dimension
%     niter : Max. number of iterations
%     tol : tolerance for the stop through the norm of the residual
%     LUv  : values of LU decomposition matrix (for SSOR or ILU(0)) in
%              CSR_packed storage, or values of diagonal of A if
%              (abs)Diagonal preconditioning.
%     pc  : preconditioning type
%       - 0 : No preconditioning
%       - 1 : Diagonal preconditioning
%       - 2 : Absolute diagonal preconditioning
%       - 3 : Symmetric SOR (Symmetric Gauss-Seidel -> SSOR, w=1)
%       - 4 : ILU(0)
%
%      Sergio A. Castiblanco B.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = max(Ac);
nb = length(b);
if n~=nb
    t = -1;
    return
end

% NO PRECONDITIONING
if pc==0

for t = 1:niter

    %Residual
    r = b - csr_matvec(Av,Ac,Ar,x);
    
    %Arnoldi Method with Householder
    [WH,pbet] = csr_arnoldi_householder(Av,Ac,Ar,r,m,LUv,0);
    
    %Elimination process for Hessemberg matrix H through Givens Rotations
    g = zeros(m+1,1);
    g(1) = WH(1,1);
    for i=1:m
        d = sqrt(WH(i,i+1)^2 + WH(i+1,i+1)^2);
        s = WH(i+1,i+1)/d;
        c = WH(i,i+1)/d;
        WH(i,i+1) = WH(i,i+1)*c + WH(i+1,i+1)*s;
        WH(i+1,i+1) = 0;
        for j=i+1:m
            %Updating Hessemberg Matrix
            hij = WH(i,j+1);
            hi1j = WH(i+1,j+1);
            WH(i,j+1) = c*hij + s*hi1j;
            WH(i+1,j+1) = -s*hij + c*hi1j;
        end
        %Updating RHS vector g
        g(i+1) = -s*g(i);
        g(i) = c*g(i);
    end
    
    %Backward sustitution
    for i = m:-1:1
        g(i) = ( g(i)-WH(i,i+2:m+1)*g(i+1:m) )/WH(i,i+1);
    end
    %here y = g(1:m)
    
    %Computing z for xm = x0+z
    z = zeros(n,1);
    for k = m:-1:1
        z(k) = z(k) + g(k);
        q = z;
        bqTw = pbet(k)*(q(k:n)'*WH(k+1:n+1,k));
        for i = k:n
            z(i) = q(i) - WH(i+1,k)*bqTw;
        end
    end
    
    %updating the solution
    x = x + z;
    
    %tolerance over residual
    res = abs(g(m+1));
    if res<tol
        break
    end
end

%DIAGONAL OR ABSOLUTE DIAGONAL PRECONDITIONING
elseif pc==1 || pc==2

% D = csr_diaga(Av,Ac,Ar);
for t = 1:niter

    %Residual
    r = b - csr_matvec(Av,Ac,Ar,x);
    
    %Arnoldi Method with Householder
    [WH,pbet] = csr_arnoldi_householder(Av,Ac,Ar,r,m,LUv,1);
    
    %Elimination process for Hessemberg matrix H through Givens Rotations
    g = zeros(m+1,1);
    g(1) = WH(1,1);
    for i=1:m
        d = sqrt(WH(i,i+1)^2 + WH(i+1,i+1)^2);
        s = WH(i+1,i+1)/d;
        c = WH(i,i+1)/d;
        WH(i,i+1) = WH(i,i+1)*c + WH(i+1,i+1)*s;
        WH(i+1,i+1) = 0;
        for j=i+1:m
            %Updating Hessemberg Matrix
            hij = WH(i,j+1);
            hi1j = WH(i+1,j+1);
            WH(i,j+1) = c*hij + s*hi1j;
            WH(i+1,j+1) = -s*hij + c*hi1j;
        end
        %Updating RHS vector g
        g(i+1) = -s*g(i);
        g(i) = c*g(i);
    end
    
    %Backward sustitution
    for i = m:-1:1
        g(i) = ( g(i)-WH(i,i+2:m+1)*g(i+1:m) )/WH(i,i+1);
    end
    %here y = g(1:m)
    
    %Computing z for xm = x0+z
    z = zeros(n,1);
    for k = m:-1:1
        z(k) = z(k) + g(k);
        q = z;
        bqTw = pbet(k)*(q(k:n)'*WH(k+1:n+1,k));
        for i = k:n
            z(i) = q(i) - WH(i+1,k)*bqTw;
        end
    end
    z = z./LUv;       %preconditioning change
    
    %updating the solution
    x = x + z;
    
    %tolerance over residual
    res = abs(g(m+1));
    if res<tol
        break
    end
end

%SSOR OR ILU(0) PRECONDITIONING
elseif pc==3 || pc==4

for t = 1:niter

    %Residual
    r = b - csr_matvec(Av,Ac,Ar,x);
    
    %Arnoldi Method with Householder
    [WH,pbet] = csr_arnoldi_householder(Av,Ac,Ar,r,m,LUv,3);
    
    %Elimination process for Hessemberg matrix H through Givens Rotations
    g = zeros(m+1,1);
    g(1) = WH(1,1);
    for i=1:m
        d = sqrt(WH(i,i+1)^2 + WH(i+1,i+1)^2);
        s = WH(i+1,i+1)/d;
        c = WH(i,i+1)/d;
        WH(i,i+1) = WH(i,i+1)*c + WH(i+1,i+1)*s;
        WH(i+1,i+1) = 0;
        for j=i+1:m
            %Updating Hessemberg Matrix
            hij = WH(i,j+1);
            hi1j = WH(i+1,j+1);
            WH(i,j+1) = c*hij + s*hi1j;
            WH(i+1,j+1) = -s*hij + c*hi1j;
        end
        %Updating RHS vector g
        g(i+1) = -s*g(i);
        g(i) = c*g(i);
    end
    
    %Backward sustitution
    for i = m:-1:1
        g(i) = ( g(i)-WH(i,i+2:m+1)*g(i+1:m) )/WH(i,i+1);
    end
    %here y = g(1:m)
    
    %Computing z for xm = x0+z
    z = zeros(n,1);
    for k = m:-1:1
        z(k) = z(k) + g(k);
        q = z;
        bqTw = pbet(k)*(q(k:n)'*WH(k+1:n+1,k));
        for i = k:n
            z(i) = q(i) - WH(i+1,k)*bqTw;
        end
    end
    z = csr_solpacklu(LUv, Ac, Ar, z);       %preconditioning change
    
    %updating the solution
    x = x + z;
    
    %tolerance over residual
    res = abs(g(m+1));
    if res<tol
        break
    end
end

end

end