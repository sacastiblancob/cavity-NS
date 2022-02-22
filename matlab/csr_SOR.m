function [x,t] = csr_SOR(Av,Ar,Ac,Pv,Pr,Pc,Qv,Qr,Qc,b,x,niter,tol)
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
%     Av, Ac, Ar : Components of matrix A in CSR stotage
%     Pv, Pc, Pr : Components of matrix P in CSR stotage
%     Qv, Qc, Qr : Components of matrix Q in CSR stotage
%     b : right hand side vector
%     x : first guess for the solution of the system
%     niter : max. number of iterations for solve the system
%     tol : tolerance for the residual, if the method reaches the tol value
%           then its finished
%
%%%%%

n = length(Ar)-1;

% loop
for t = 1:niter
    x = csr_matvec(Qv,Qc,Qr,x) + b;
    for i=1:n
        j = Pr(i):Pr(i+1)-Pr(1);
        if length(j)==1
            x(i) = x(i)/Pv(j(length(j)));
        else
            x(i) = (x(i) - Pv(j(1:length(j)-1))*x(Pc(j(1:length(j)-1))))/Pv(j(length(j)));
        end
    end
    if norm(b - csr_matvec(Av,Ac,Ar,x))<tol
        return
    end
end

end