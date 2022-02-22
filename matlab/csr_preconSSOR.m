function [PQv] = csr_preconSSOR(Av,Ac,Ar,w)
%
% This function takes the matrix A in CSR storage, and prepares the
% matrices for solve the system Ax=b with SSOR preconditioner
%
% P and Q:
% The SSOR preconditioner is defined with the next matrix:
% 
%               Mssor = (D - w*E) D^(-1) (D - wF)
%
% In which the system can be solved throug the matrices:
%
%               P = (D- w*E)D^(-1) = (I - w*E*D^(-1))
%               Q = (D - w*F)
%
% Then, for finding z = Mssor^(-1)*x
%
%               solve Py = x;  -> Forward substitution
%               solve Qz = y;  -> Backward substitution
%
% For the particular case when w = 1, we are in Symmetric Gauss-Seidel
% preconditioner.
% Note that diagonal of matrix P = 1. Then this decomposition on P and Q
% can be stored in packed format,i.e., use the same structure as A, but
% storing in the lower part -w*E*D^(-1), and in the upper part Q.
%
% Entries:
%     Av, Ac, Ar : Components of matrix A in CSR stotage
%     w : SOR Over-Relaxation coefficient
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(Ar)-1;
nz = length(Av);
D = zeros(n,1);

for j=1:n
    for i=Ar(j):Ar(j+1)-1
        if Ac(i) == j
            D(j) = Av(i);
        end
    end
end

PQv = zeros(1,nz);

pqp = 1;
for j=1:n
    for i=Ar(j):Ar(j+1)-1
        if Ac(i) > j
            PQv(pqp) = w*Av(i);
        elseif Ac(i) < j
            PQv(pqp) = w*Av(i)/D(Ac(i));
        else
            PQv(pqp) = Av(i);
        end
        pqp = pqp + 1;
    end
end

end