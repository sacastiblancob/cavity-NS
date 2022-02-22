function [LUv] = csr_preconILU0(Av,Ac,Ar)
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
% Original algorithm: (Mittal & Al-Kurdi, 2003)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LUv = Av;

%Incomplete LU sparse CSR (Mittal Al-Kurdi, 2003)
n = length(Ar)-1;
dia = zeros(n,1);
for i=1:n
    for j=Ar(i):Ar(i+1)-1
        if Ac(j) == i
            dia(i) = j;
        end
    end
end
point = zeros(n,1);

for i = 2:n
    for v = Ar(i)+1:Ar(i+1)-1
        point(Ac(v))=v;
    end
    for v = Ar(i):dia(i)-1
        j = Ac(v);
        LUv(v) = LUv(v)/LUv(dia(j));    
        for w = dia(j)+1:Ar(j+1)-1
            k = point(Ac(w));
            if k>0
                LUv(k) = LUv(k) - LUv(v)*LUv(w);
%             end
            %else
            %   LUv(dia(i)) = LUv(dia(i)) - LUv(v)*LUv(w);
            end
        end
    end
    for v=Ar(i)+1:Ar(i+1)-1
        point(Ac(v))=0;
    end
end

end