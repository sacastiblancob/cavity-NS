function [Pv] = csc_diaga(Av,Ar,Ac)
%
% This function takes the matrix A in CSC storage, and gets the diagonal of
% A.
%
% Entries are the values(Av), rows(Ar) and columns(Ac) in CSC storage of A
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%
%

m = length(Ac)-1;
Pv = zeros(m,1);
for j=1:m
    for i=Ac(j):Ac(j+1)-1
        if Ar(i) == j
            Pv(j) = Av(i);
        end
    end
end

end