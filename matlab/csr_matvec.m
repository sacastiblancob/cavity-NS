function c = csr_matvec(Av,Ac,Ar,b)
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
m = length(Ar)-1;
n = m;
if m~=length(b)
    disp('ERROR!!! Dimensions does not agree')
    return
end
c = zeros(n,1);

for i=1:n
    k1 = Ar(i);
    k2 = Ar(i+1)-1;
    c(i) = Av(k1:k2)*b(Ac(k1:k2));
end

end