function [V,H] = csr_arnoldi_MGS(Av,Ac,Ar,v,m,epsro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes de Arnoldi's Method with Modified Gram-Schmidt
%
% Entries:
%   - Av,Ac,Ar : matrix A of dimensions n*n in CSR storage
%   - v : vector for compute Krylov subspace = span{v,Av,(A^2)v,...,(A^m-1)v}
%   - m : dimension of Krylov subspace
%   - epsro : epsilon threshold for reorthogonalization
%
% Output
%   - V (n*n) : orthonormal base of Krylov Subspace
%   - H (m*m) : A-similar Hessemberg Matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = v./norm(v);
n = length(Ac) - 1;

H = zeros(m,m);
V = zeros(n,m);
V(:,1) = v;
for j=1:m
    w = csr_matvec(Av,Ac,Ar,V(:,j));
    w0 = w;     %for comparing when reorthogonalization
    for i=1:j
        H(i,j) = w'*V(:,i);
        w = w - H(i,j)*V(:,i);
    end
    if abs(norm(w) - norm(w0)) < epsro %reorthogonalization if needed (almost never needed)
        w = csr_matvec(Av,Ac,Ar,V(:,j));
        w = -w;
        for i=1:j
            H(i,j) = w'*V(:,i);
            w = w - H(i,j)*V(:,i);
        end
    end
    if j~=m
        H(j+1,j) = norm(w);
        if H(j+1,j) <= 2*eps(1)
            break
        end   
        V(:,j+1) = w/H(j+1,j);
    end
end

end