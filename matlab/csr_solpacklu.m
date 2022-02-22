function x = csr_solpacklu(LUv, LUc, LUr, b)
% This function compute the forward and backward sustitution for the system
% LUx = b, where LU should be entered in pack_csr notation, that's it:
% packLU in Compressed Sparse Row notation.
% Enries are:
%       LUv = vector with values of packLU matrix in CSR format
%       LUc = vector with column indexes of packLU matrix in CSR format
%       LUr = vector with row pointers of packLU matrix in CSR format

x = b;
n = length(LUr)-1;

%forward sustitution
for i=1:n
    for j = LUr(i):LUr(i+1)-1
        if i>LUc(j)
            x(i) = x(i) - LUv(j)*x(LUc(j));
        end
    end
end

%backward sustitution
for i=n:-1:1
    for j = LUr(i+1)-1:-1:LUr(i)
        if i<LUc(j)
            x(i) = x(i) - LUv(j)*x(LUc(j));
        elseif i==LUc(j)
            x(i) = x(i)/LUv(j);
        end
    end
end

end