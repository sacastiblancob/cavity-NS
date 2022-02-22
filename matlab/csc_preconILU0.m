function [LUv] = csc_preconILU0(Av,Ar,Ac)
%
% This function takes the matrix A in CSC storage, and prepares the
% matrices for solve the system Ax=b with ILU0 preconditioner
%
% L and U:
% The ILU0 preconditioner is defined with the next matrix:
% 
%               Milu0 = L*U
%
% Where L and U are the result of an Incomplete LU factorization, with the
% same non-zero structure as matrix A. Since L is a unitary lower
% triangular matrix. The ILU0 can be stored in packed format, and, as it
% has the same structure as A, it is only needed to compute the vector of
% values LUv. Additionaly
%
% In which the system can be solved in the next way:
%
%               solve Ly = x;  -> Forward substitution
%               solve Uz = y;  -> Backward substitution
%
% Entries:
%     Av, Ar, Ac : Components of matrix A in CSC stotage
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = length(Ac)-1;
LUv = Av;

%The algorithm for Incomplete LU
for j=1:m-1
    for i=Ac(j):Ac(j+1)-1
        if j == Ar(i)
            d = LUv(i);
        end
        if j < Ar(i)
            LUv(i) = LUv(i)/d;
            for j2 = j+1:m
                a = 0;
                for i2 = Ac(j2):Ac(j2+1)-1
                    if Ar(i2)==j
                       a = LUv(i2);
                    end
                    if Ar(i)==Ar(i2)
                        LUv(i2) = LUv(i2) - LUv(i)*a;
                    end
                end
            end            
        end
    end
end

end