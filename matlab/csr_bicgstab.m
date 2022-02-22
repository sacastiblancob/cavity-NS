function [x,t,res] = csr_bicgstab(Av,Ac,Ar,b,x,r0a,niter,tol,LUv,pc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function solves the system Ax=b with the Biconjugate Gradient
% Stabilized Method (Van der Vorst algorithm).
%
% Entries:
%     Av,Ac,Ar : Matrix A in CSR storage
%     b : Right hand side vector
%     x : First guest for the solution
%     r0a : Arbitrary r0 for computing conjugate directions
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
% Preconditioning is based on Analyisis and Practical use of Flexible
% BICGSTAB (Jie Chen, Louis Curfman, Hong Zhang)
%
%      Sergio A. Castiblanco B.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%No preconditioning
if pc==0

r = b - csr_matvec(Av,Ac,Ar,x);
p = r;
for t=1:niter
    Ap = csr_matvec(Av,Ac,Ar,p);
    a = (r'*r0a)/(Ap'*r0a);
    s = r - a*Ap;
    As = csr_matvec(Av,Ac,Ar,s);
    w = (As'*s)/(As'*As);
    x = x + a*p + w*s;
    b = (a/w)*((s - w*As)'*r0a)/(r'*r0a);
    r = s - w*As;
    if norm(r)<tol
        break
    end
    p = r + b*(p - w*Ap);
end
res = norm(r);

%Diagonal or absolute diagonal preconditioning
elseif pc==1 || pc==2
    
r = b - csr_matvec(Av,Ac,Ar,x);
p = r;
% D = csr_diaga(Av,Ac,Ar);
for t=1:niter
    pg = p./LUv;
    Ap = csr_matvec(Av,Ac,Ar,pg);
    a = (r'*r0a)/(Ap'*r0a);
    s = r - a*Ap;
    sg = s./LUv;
    As = csr_matvec(Av,Ac,Ar,sg);
    w = (As'*s)/(As'*As);
    x = x + a*pg + w*sg;
    b = (a/w)*((s - w*As)'*r0a)/(r'*r0a);
    r = s - w*As;
    if norm(r)<tol
        break
    end
    p = r + b*(p - w*Ap);
end
res = norm(r);  

%SSOR or ILU(0) preconditioning
elseif pc==3 || pc==4

r = b - csr_matvec(Av,Ac,Ar,x);
p = r;
for t=1:niter
    pg = csr_solpacklu(LUv, Ac, Ar, p);
    Ap = csr_matvec(Av,Ac,Ar,pg);
    a = (r'*r0a)/(Ap'*r0a);
    s = r - a*Ap;
    sg = csr_solpacklu(LUv, Ac, Ar, s);
    As = csr_matvec(Av,Ac,Ar,sg);
    w = (As'*s)/(As'*As);
    x = x + a*pg + w*sg;
    b = (a/w)*((s - w*As)'*r0a)/(r'*r0a);
    r = s - w*As;
    if norm(r)<tol
        break
    end
    p = r + b*(p - w*Ap);
end
res = norm(r);  

end

end