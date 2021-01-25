function [x,t] = csc_SOR(Av,Ar,Ac,Pv,Pr,Pc,Qv,Qr,Qc,b,x,niter,tol)
%Solves the system Ax = b with the SOR iterative method (Symmetric Over 
% Relaxation).
%
% P and Q:
% Classic iterative methods use the next recurrence equation for solve the 
% system:
% 
%               P * x(t+1) = Q * x(t) + b
%
%       For SOR:
%            P = ((1/w)*D + L)
%            Q = ((1/w - 1)*D - U); where: L = lower part of A
%                                          U = upper part of A
%                                          D = diagonal of A
%
% For the particular case when w = 1, we are in Gauss-Seidel method
%
% Entries:
%     Av, Ar, Ac : Components of matrix A in CSC stotage
%     Pv : Vector with the values of D (diagonal of A)
%     Qv, Qr, Qc : Components of matrix Q in CSC stotage
%     b : right hand side vector
%     x : first guess for the solution of the system
%     niter : max. number of iterations for solve the system
%     tol : tolerance for the residual, if the method reaches the tol value
%           then its finished
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%      Pontificia Universidad Javeriana - Bogotá
%

m = length(Ac)-1;

% loop
for t = 1:niter
    x = csc_matvec(Qv,Qr,Qc,x) + b;
    for j=1:m
        for i = Pc(j):Pc(j+1)-Pc(1)
            if j==Pr(i)
                x(j) = x(j)/Pv(i);
            else
                x(Pr(i)) = x(Pr(i)) - Pv(i)*x(j);
            end
        end
    end
    if norm(b - csc_matvec(Av,Ar,Ac,x))<tol
        return
    end
end

end