function [Cc,Cr] = csr_symmatmul(Ac,Ar,Bc,Br)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the symbolic multiplication of matrices stored in
% Compresed Sparse Row (CSR). Vectors Ac,Bc,Cc are the column indices and
% the vectors Cr,Ar,Br are the row pointers respectively.
%
%    C = A*B
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of rows of the first matrix
np = length(Ar)-1;

%Number of columns of the first matrix and number of rows of second matrix
nq1 = max(Ac);
nq2 = length(Br)-1;
if nq1~=nq2
    Cc = NaN;
    Cr = NaN;
    return
end
nq = nq1;

%number of columns of the second matrix
nr = max(Bc);

%allocating Cr
Cr = zeros(1,np+1);

%working array for multiple switch
ix = zeros(1,nr);

%row pointer counter
ip = 1;

%computing length of Cc
for i=1:np
    Cr(i) = ip;
    iaa = Ar(i);
    iab = Ar(i+1)-1;
    for jp = iaa:iab
        j = Ac(jp);
        iba = Br(j);
        ibb = Br(j+1)-1;
        for kp = iba:ibb
            k = Bc(kp);
            if ix(k)~=i
%                 Cc(ip) = k;
                ip = ip + 1;
                ix(k) = i;
            end
        end
    end
end
Cr(np+1) = ip;

%Allocating Cc
Cc = zeros(1,ip-1);

%%%Filling Cc
%row pointer counter
ip = 1;

%multiple switch vector
ix = zeros(1,nr);

for i=1:np
%     Cr(i) = ip;
    iaa = Ar(i);
    iab = Ar(i+1)-1;
    for jp = iaa:iab
        j = Ac(jp);
        iba = Br(j);
        ibb = Br(j+1)-1;
        for kp = iba:ibb
            k = Bc(kp);
            if ix(k)~=i
                Cc(ip) = k;
                ip = ip + 1;
                ix(k) = i;
            end
        end
    end
end

end