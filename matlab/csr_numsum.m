function [Cv,Cc,Cr] = csr_numsum(Av,Ac,Ar,Bv,Bc,Br,Cc,Cr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes numerical sum of matrices stored in Compresed
% Sparse Row (CSR). Vectors Av,Ac,Ar are the matrix A stored in CSR storage
% Bv,Bc,Br are the matrix B stored in CSR storage. Cc,Cr are the structure
% of matrix C computed previously with csr_symsum
%
%    C = A + B
%
% * if a cero appears due to A(i,j) + B(i,j) = 0; it will be preserved as
%   a non-zero, even when strictly speaking it is not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(Ar)-1;
m = max(max(Ac),max(Bc));

%Allocating Cv
Cv = zeros(1,Cr(n+1)-1);

%Working array
x = zeros(1,m);

for i=1:n
    ica = Cr(i);
    icb = Cr(i+1)-1;
    for ip=ica:icb
        x(Cc(ip))=0;
    end
    iaa = Ar(i);
    iab = Ar(i+1)-1;
    for ip = iaa:iab
        x(Ac(ip))=Av(ip);
    end
    iba = Br(i);
    ibb = Br(i+1)-1;
    for ip = iba:ibb
        j = Bc(ip);
        x(j) = x(j) + Bv(ip);
    end
    for ip=ica:icb
        Cv(ip) = x(Cc(ip));
    end    
end

end