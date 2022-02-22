function [x,t] = csc_DPCG(Av,Ar,Ac,b,x,niter,tol)
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
m = length(Ac)-1;

P = zeros(m,1);
for j=1:m
    for i=Ac(j):Ac(j+1)-1
        if Ar(i) == j
            P(j) = Av(i);
        end
    end
end

ro = b - csc_matvec(Av,Ar,Ac,x);

zo = ro./P;

d = zo;
for t = 1:niter
    Ad = csc_matvec(Av,Ar,Ac,d);
    alf = (ro'*zo)/(d'*Ad);
    x = x + alf*d;
    if norm(b - csc_matvec(Av,Ar,Ac,x))<tol
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