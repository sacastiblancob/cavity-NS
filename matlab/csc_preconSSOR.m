function [PQv] = csc_preconSSOR(Av,Ar,Ac,w)
%
% This function takes the matrix A in CSC storage, and prepares the
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
%     Av, Ar, Ac : Components of matrix A in CSC stotage
%     w : SOR Over-Relaxation coefficient
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = length(Ac)-1;
nz = length(Av); 
D = zeros(m,1);

for j=1:m
    for i=Ac(j):Ac(j+1)-1
        if Ar(i) == j
            D(j) = Av(i);
        end
    end
end

PQv = zeros(nz,1);

pqp=1;
for j=1:m
    for i=Ac(j):Ac(j+1)-1
        if Ar(i) < j
            PQv(pqp) = w*Av(i);
        elseif Ar(i) > j
            PQv(pqp) = w*Av(i)/D(j);
        else
            PQv(pqp) = Av(i);
        end
        pqp = pqp + 1;
    end
end

end