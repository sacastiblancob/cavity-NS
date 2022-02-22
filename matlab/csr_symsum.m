function [Cc,Cr] = csr_symsum(Ac,Ar,Bc,Br)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the symbolic sum of matrices stored in Compresed
% Sparse Row (CSR). Vectors Ac,Bc,Cc are the column indices and the vectors
% Cr,Ar,Br are the row pointers respectively.
%
%    C = A + B
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

an = length(Ar)-1;
bn = length(Br)-1;
if an~=bn
    disp('Dimensions does not agree!!!, errors might be expected!!!')
    Cc = NaN;
    Cr = NaN;
    return
end
n = an;
m = max(max(Ac),max(Bc));

%row pointer counter
ip = 1;

%alocating Cr
Cr = zeros(1,n+1);

%multiple switch vector
ix = zeros(1,m);

for i=1:n
   Cr(i)=ip;
   iaa = Ar(i);
   iab = Ar(i+1)-1;
   for jp = iaa:iab
       j = Ac(jp);
%        Cc(ip) = j;
       ip = ip + 1;
       ix(j) = i;      
   end
   iba = Br(i);
   ibb = Br(i+1)-1;
   for jp = iba:ibb
       j = Bc(jp);
       if (ix(j)==i)
           continue
       else
%            Cc(ip) = j;
           ip = ip + 1;
       end
   end
end
Cr(n+1) = ip;

%Allocating Cc
Cc = zeros(1,ip-1);

%%%Filling Cc
%row pointer counter
ip = 1;

%multiple switch vector
ix = zeros(1,m);

for i=1:n
%    Cr(i)=ip;
   iaa = Ar(i);
   iab = Ar(i+1)-1;
   for jp = iaa:iab
       j = Ac(jp);
       Cc(ip) = j;
       ip = ip + 1;
       ix(j) = i;      
   end
   iba = Br(i);
   ibb = Br(i+1)-1;
   for jp = iba:ibb
       j = Bc(jp);
       if (ix(j)==i)
           continue
       else
           Cc(ip) = j;
           ip = ip + 1;
       end
   end
end

end