function [Tv,Tc,Tr] = css_tperm(Av,Ac,Ar,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes the matrix A in CSR storage, and computes the
% transpose of the permutation given on vector w such as:
%
%      (PA)'  ; P = permutation matrix of p
%
% If you want to apply the permutation over rows and columns, you need to
% use this algorithm twice, obtaining
%
%      (PAQ)' ; P = permutation matrix of p; Q = permutation matrix of q
%               where p is the vector for permutations of row
%                     q is the vector for permutations of columns
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of columns and non-zeros of A
n = length(Ar)-1;
nz = Ar(n+1)-1;

%Number of columns of T
mt = max(Ac);

%Allocating Tv, Tr and Tc
Tv = zeros(1,nz);
Tc = zeros(1,nz);
Tr = zeros(1,n+1);

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
    iaa = Ar(p(i));
    iab = Ar(p(i)+1)-1;
    if iab<iaa
        break
    end
    for jp=iaa:iab
        j=Ac(jp)+1;
        k=Tr(j);
%         Tc(k)=p(i); typo in the book (Sparse Matrix Technology, chap. 7.5)
        Tc(k)=i;
        Tv(k)=Av(jp);
        Tr(j)=k+1;
    end
end    

end