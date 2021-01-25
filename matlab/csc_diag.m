function [Av,Ar,Ac] = csc_diag(h2,d)
%
%This function make a matrix in CSC storage by diagonals, the entries are:
%      h2 --> matrix whose columns are gonna be the diagonals of A
%      d ---> row-vector with the numbers of the diagonals in h2 cols
%             0 for the main diagonal, >0 for lower part and >0 for upper
%             none of the diagonals are mandatory
%
%      warning: if you have zeros in h2 you should remove them with another
%          function.
%      NOTE: This function computes a square matrix with the given
%      diagonals.
%
%      Example:
%                    h2              d
%               | 6  1  6 |      [-1 0 1]
%               | 7  2  7 |
%               | 8  3  8 |
%               | 9  4  9 |
%               | 10 5 10 |
%      Returns:
%                       A
%              | 1  6  0  0  0 |
%              | 6  2  7  0  0 |
%              | 0  7  3  8  0 |
%              | 0  0  8  4  9 |
%              | 0  0  0  9  5 |
%
%      Note that 10's of h2(5,1) and h2(5,3) have dissapeared, those could
%      be zero anyway
%
%      Of course A is stored in CSC, so the function returns
%      Av, Ar, Ac (values, rows, columns)
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%      Pontificia Universidad Javeriana - Bogotá
%

siz = size(h2);
if siz(2)~=length(d)
    disp('Dimensions of h2 and d does not agree!!')
    return
end

dp = d(abs(d) == min(abs(d)));
dp = dp(1);
m = length(h2(:,d==dp))+abs(dp);
nz=0;
for i=1:length(d)
    %nz = nz+length(h2)-abs(d(i))+(m-length(h2));
    nz = nz - abs(d(i)) + m;
end
Av = zeros(nz,1);
Ar = zeros(nz,1); 
Ac = zeros(m+1,1);
Ac(1) = 1;
k = 0;
c = -m;
cv = d*0;
p = 1;
for i=1:m
    Ac(i+1) = Ac(i) + sum((d<=k).*(d>c));
    Ar(Ac(i):Ac(i+1)-1) = sort(abs(d(logical((d<=k).*(d>c)))-1-k));
    cv = cv+(d<=k).*(d>c);
    hcols = h2(:,logical((d<=k).*(d>c)));
    hrows = flip(cv(logical((d<=k).*(d>c))));
    row = 1;
    col = length(hrows);
    for j=Ac(i):Ac(i+1)-1
        Av(p) = hcols(hrows(row),col);
        row = row + 1;
        col = col - 1;
        p = p + 1;
    end
    k=k+1;
    c=c+1;
end

end