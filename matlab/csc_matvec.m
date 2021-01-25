function c = csc_matvec(Av,Ar,Ac,b)
%
% This function computes the multiplication matrix-vector with the matrix
% A stored in CSC
%
%    Ab=c
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%      Pontificia Universidad Javeriana - Bogotá
%

%n = max(Ar);
m = length(Ac)-1;
n = m;
% if m~=length(b)
%     disp('ERROR!!! Dimensions does not agree')
%     return
% end
c = zeros(n,1);

for j=1:m
    for i=Ac(j):Ac(j+1)-1
       r = Ar(i);
       c(r) = c(r) + Av(i)*b(j);
    end
end

end