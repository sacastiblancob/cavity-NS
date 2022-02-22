function [Tc,Tr] = css_symtrans(Ac,Ar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes the matrix A in CS storage, and computes its
% symbolic transpose such as:
%
%      T = A'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of columns and non-zeros of A
n = length(Ar)-1;
nz = Ar(n+1)-1;

%Number of columns of T
mt = max(Ac);

%Allocating Tr and Tc
Tc = zeros(1,nz);
Tr = zeros(1,mt+1);

mh = mt+1;

for i=1:nz
    j = Ac(i)+2;
    if(j<=mh)
        Tr(j)=Tr(j)+1;
    end
end

Tr(1) = 1;
Tr(2) = 1;

for i=3:mh
    Tr(i)=Tr(i)+Tr(i-1);
end

for i=1:n
    iaa = Ar(i);
    iab = Ar(i+1)-1;
    if iab<iaa
        break
    end
    for jp=iaa:iab
        j=Ac(jp)+1;
        k=Tr(j);
        Tc(k)=i;
        Tr(j)=k+1;
    end
end

end