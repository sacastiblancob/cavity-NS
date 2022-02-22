function [x,t,res] = csc_bicgstab(Av,Ar,Ac,b,x,r0a,niter,tol,LUv,pc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function solves the system Ax=b with the Biconjugate Gradient
% Stabilized Method (Van der Vorst algorithm).
%
% Entries:
%     Av,Ar,Ac : Matrix A in CSC storage
%     b : Right hand side vector
%     x : First guest for the solution
%     r0a : Arbitrary r0 for computing conjugate directions
%     niter : Max. number of iterations
%     tol : tolerance for the stop through the norm of the residual
%     LUv  : values of LU decomposition matrix (for SSOR or ILU(0)) in
%              CSC_packed storage, if not SSOR or ILU(0), it must be have
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

r = b - csc_matvec(Av,Ar,Ac,x);
p = r;
for t=1:niter
    Ap = csc_matvec(Av,Ar,Ac,p);
    a = (r'*r0a)/(Ap'*r0a);
    s = r - a*Ap;
    As = csc_matvec(Av,Ar,Ac,s);
    w = (As'*s)/(As'*As);
    x = x + a*p + w*s;
    bet = (a/w)*((s - w*As)'*r0a)/(r'*r0a);
    r = s - w*As;
    if norm(r)<tol
        break
    end
    p = r + bet*(p - w*Ap);
end
res = norm(r);

%Diagonal or Absolute Diagonal preconditioning
elseif pc==1 || pc==2
    
r = b - csc_matvec(Av,Ar,Ac,x);
p = r;
% D = csc_diaga(Av,Ar,Ac);
for t=1:niter
    pg = p./LUv;
    Ap = csc_matvec(Av,Ar,Ac,pg);
    a = (r'*r0a)/(Ap'*r0a);
    s = r - a*Ap;
    sg = s./LUv;
    As = csc_matvec(Av,Ar,Ac,sg);
    w = (As'*s)/(As'*As);
    x = x + a*pg + w*sg;
    bet = (a/w)*((s - w*As)'*r0a)/(r'*r0a);
    r = s - w*As;
    if norm(r)<tol
        break
    end
    p = r + bet*(p - w*Ap);
end
res = norm(r);  

%SSOR or ILU(0) preconditioning
elseif pc==3 || pc==4

r = b - csc_matvec(Av,Ar,Ac,x);
p = r;
for t=1:niter
    pg = csc_solpacklu(LUv, Ar, Ac, p);
    Ap = csc_matvec(Av,Ar,Ac,pg);
    a = (r'*r0a)/(Ap'*r0a);
    s = r - a*Ap;
    sg = csc_solpacklu(LUv, Ar, Ac, s);
    As = csc_matvec(Av,Ar,Ac,sg);
    w = (As'*s)/(As'*As);
    x = x + a*pg + w*sg;
    bet = (a/w)*((s - w*As)'*r0a)/(r'*r0a);
    r = s - w*As;
    if norm(r)<tol
        break
    end
    p = r + bet*(p - w*Ap);
end
res = norm(r);  

end

end