function [v,c,r] = full2csr(A)
% This function stores the matrix A in CSC (Compressed Sparse Column),
% could be more efficient using some special matlab tools, but is made it
% with the purpouse to do it as well in fortran.
% Returns an structure with next elements:
%   v = row-vector with values of non-zero elements
%   r = row-vector with col-index
%   c = row-vector with row-index in v

    siz = size(A);
    n = siz(1);
    m = siz(2);
    nz=0;

    for j=1:m
        for i=1:n
            if A(i,j)~=0
                nz=nz+1;
            end
        end
    end

    v = zeros(1,nz);      %Vector to store values
    c = zeros(1,nz);      %Vector to store the row indexes
    r = zeros(1,n+1);     %Vector to store where column begins
    r(1) = 1;

    k=1;
    for i=1:n
        for j=1:m
            if A(i,j)~=0
                v(k) = A(i,j);
                c(k) = j;
                k = k+1;
            end
        end
        r(i+1) = k;
    end

end