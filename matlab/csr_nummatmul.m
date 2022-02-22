function [Cv,Cc,Cr] = csr_nummatmul(Av,Ac,Ar,Bv,Bc,Br,Cc,Cr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes numerical multiplication of matrices stored in
% Compresed Sparse Row (CSR). Vectors Av,Ac,Ar are the matrix A stored in
% CSR storage Bv,Bc,Br are the matrix B stored in CSR storage. Cc,Cr are
% the structure of matrix C computed previously with csr_symsum
%
%    C = A*B
%
% * if a cero appears due to A(i,*)*B(*,j) = 0; it will be preserved as
%   a non-zero, even when strictly speaking it is not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of rows of A (i.e. number of rows of C)
np = length(Ar)-1;

%number of columns of the second matrix
nr = max(Bc);

%working array
x = zeros(1,nr);

%allocating Cv
Cv = zeros(1,length(Cc));

for i=1:np
    ica = Cr(i);
    icb = Cr(i+1)-1;
    for j=ica:icb
        x(Cc(j)) = 0;
    end
    iaa = Ar(i);
    iab = Ar(i+1)-1;
    for jp = iaa:iab
        j = Ac(jp);
        a = Av(jp);
        iba = Br(j);
        ibb = Br(j+1)-1;
        for kp = iba:ibb
            k = Bc(kp);
            x(k) = x(k) + a*Bv(kp);
        end
    end
    for j=ica:icb
        Cv(j) = x(Cc(j));
    end
end

end