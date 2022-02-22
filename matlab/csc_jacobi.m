function [x,t] = csc_jacobi(Av,Ar,Ac,Pv,Qv,Qr,Qc,b,x,niter,tol)
%Solves the system Ax = b with the Jacobi iterative method.
%
% P and Q:
% Classic iterative methods use the next recurrence equation for solve the 
% system:
% 
%               P * x(t+1) = Q * x(t) + b
%
%       For jacobi:
%            P = D = diagonal of matrix A
%            Q = -(L + U); where: L = lower part of A
%                                 U = upper part of A
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

%loop
for t=1:niter
    d = csc_matvec(Qv,Qr,Qc,x) + b;
    x = d./Pv;
    if norm(b - csc_matvec(Av,Ar,Ac,x))<tol
        return
    end
end

end