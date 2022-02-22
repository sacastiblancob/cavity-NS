function [x,t] = csr_DPCG(Av,Ac,Ar,b,x,niter,tol)
%
%This function solves the system Ax=b with the Diagonal Preconditioned
%Conjugate Gradient method
%
% Entries:
%     Av,Ar,Ac : Matrix A in CSC storage
%     b : Right hand side vector
%     x : First guest for the solution
%     niter : Max. number of iterations
%     tol : tolerance for the stop through the norm of the residual
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%      Pontificia Universidad Javeriana - Bogotá
%


% Preconditioned CG CSC
n = length(Ar)-1;

P = zeros(n,1);
for i=1:n
    for j=Ar(i):Ar(i+1)-1
        if Ac(j) == i
            P(i) = Av(j);
        end
    end
end

ro = b - csr_matvec(Av,Ac,Ar,x);

zo = ro./P;

d = zo;
for t = 1:niter
    Ad = csr_matvec(Av,Ac,Ar,d);
    alf = (ro'*zo)/(d'*Ad);
    x = x + alf*d;
    if norm(b - csr_matvec(Av,Ac,Ar,x))<tol
        return
    end
    r = ro - alf*Ad;

    z = r./P;
    
    bet = (r'*z)/(ro'*zo);
    d = z + bet*d;
    ro = r;
    zo = z;
end

end