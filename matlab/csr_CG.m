function [x,t] = csr_CG(Av,Ac,Ar,b,x,niter,tol,LUv,pc)
%
%This function solves the system Ax=b with the Conjugate Gradient method
%
% Entries:
%     Av,Ar,Ac : Matrix A in CSC storage
%     b : Right hand side vector
%     x : First guest for the solution
%     niter : Max. number of iterations
%     tol : tolerance for the stop through the norm of the residual
%     LUv  : values of LU decomposition matrix (for SSOR or ILU(0)) in
%              CSR_packed storage, if not SSOR or ILU(0), it must be have
%              arbitrary values.
%     pc  : preconditioning type
%       - 0 : No preconditioning
%       - 1 : Diagonal preconditioning
%       - 2 : Absolute diagonal preconditioning
%       - 3 : Symmetric SOR (Symmetric Gauss-Seidel -> SSOR, w=1)
%       - 4 : ILU(0)
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%No preconditioning
if pc==0

r = b - csr_matvec(Av,Ac,Ar,x);
p = r;
for t = 1:niter
    rr = r'*r;
    Ap = csr_matvec(Av,Ac,Ar,p);
    a = rr/(Ap'*p);
    x = x + a*p;
    r = r - a*Ap;
    if norm(r)<tol
        break
    end
    b = (r'*r)/rr;
    p = r + b*p;
end

%Diagonal preconditioning
elseif pc==1 || pc==2

r = b - csr_matvec(Av,Ac,Ar,x);
% D = csr_diaga(Av,Ac,Ar);
z = r./LUv;
p = z;
for t=1:niter
    Ap = csr_matvec(Av,Ac,Ar,p);
    rz = r'*z;
    a = (rz)/(Ap'*p);
    x = x + a*p;
    r = r - a*Ap;
    if norm(r)<tol
        break
    end
    z = r./LUv;
    b = (r'*z)/(rz);
    p = z + b*p;
end

elseif pc==3 || pc==4

r = b - csr_matvec(Av,Ac,Ar,x);
z = csr_solpacklu(LUv, Ac, Ar, r);
p = z;
for t=1:niter
    Ap = csr_matvec(Av,Ac,Ar,p);
    rz = r'*z;
    a = (rz)/(Ap'*p);
    x = x + a*p;
    r = r - a*Ap;
    if norm(r)<tol
        break
    end
    z = csr_solpacklu(LUv, Ac, Ar, r);
    b = (r'*z)/(rz);
    p = z + b*p;
end

end

end