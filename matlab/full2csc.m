function [v,r,c] = full2csc(A)
% This function stores the matrix A in CSC (Compressed Sparse Column),
% could be more efficient using some special matlab tools, but is made it
% with the purpouse to do it as well in fortran.
% Returns an structure with next elements:
%   v = column-vector with values of non-zero elements
%   r = column-vector with row-index
%   c = column-vector with col-index in v

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

    v = zeros(nz,1);      %Vector to store values
    r = zeros(nz,1);      %Vector to store the row indexes
    c = zeros(m+1,1);     %Vector to store where column begins
    c(1) = 1;

    k=1;
    for j=1:m
        for i=1:n
            if A(i,j)~=0
                v(k) = A(i,j);
                r(k) = i;
                k = k+1;
            end
        end
        c(j+1) = k;
    end

end