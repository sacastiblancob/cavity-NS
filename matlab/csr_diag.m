function [Av,Ar,Ac] = csr_diag(h2,d)
%
%This function make a matrix in CSR storage by diagonals, the entries are:
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
%      A is stored in CSR, so the function returns
%      Av, Ac, Ar (values, columns, row_index)
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%      Pontificia Universidad Javeriana - Bogotá
%

siz = size(h2);
if siz(2)~=length(d)
    disp('Dimensions of h2 and d does not agree!!')
    return
end
%%%%%
%Index of principal diagonal
dp = d(abs(d) == min(abs(d)));
dp = flip(dp);
dp = dp(1);

%pos of principal diagonal on d
pdv = sum(d<=dp);

%Number of columns and rows
m = siz(1)+abs(dp);
n = m;

%Efective values readed from h2 (how many values per column on h2 will be keeped)
ed = m-abs(d);

%Number of nonzeros
nz = sum(ed);

%Initialization of Av, Ac, Ar
Av = zeros(1,nz);
Ac = zeros(1,nz);
Ar = ones(1,n+1);

%new h2
h2n = h2;
h2n = [h2n;zeros(n-siz(1),siz(2))];

%putting zeros atop left diagonals
for i=1:pdv-1
    h2n(:,i) = [zeros(abs(d(i)),1);h2(1:ed(i),i)];
end

% %reordered columns of h2
% roh2 = [pdv:siz(2),1:pdv-1];
% 
% %reordering h2 and d
% h2n = h2n(:,roh2);
% dn = d(roh2);

%column counter vector
cv = d+1;

%transpose of h2n and dn
siz = size(h2n);
h2n = h2n';

%filling Av, Ac, Ar
inz = 1;
for j=1:siz(1)
    for i=1:siz(2)
        if h2n(i,j)~=0
            Av(inz)=h2n(i,j);
            Ac(inz)=cv(i);
            Ar(j+1:n+1)=Ar(j+1:n+1)+1;
            inz = inz+1;
        end
    end
    cv=cv+1;
end

end