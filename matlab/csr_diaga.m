function [Pv] = csr_diaga(Av,Ac,Ar)
%
% This function takes the matrix A in CSR storage, and gets the diagonal of
% A.
%
% Entries are the values(Av), rows(Ar) and columns(Ac) in CSC storage of A
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%
%

n = length(Ar)-1;
Pv = zeros(n,1);
for i=1:n
    for j=Ar(i):Ar(i+1)-1
        if Ac(j) == i
            Pv(i) = Av(j);
        end
    end
end

end