function [v,r,c] = mat_sparse2csc(A)
% This function stores the matrix A in CSC (Compressed Sparse Column),
% could be more efficient using some special matlab tools, but is made it
% with the purpouse to do it as well in fortran.
% Returns an structure with next elements:
%   v = column-vector with values of non-zero elements
%   r = column-vector with row-index
%   c = column-vector with col-index in v

    siz = size(A);
%     n = siz(1);
    m = siz(2);

%     v = zeros(nz,1);      %Vector to store values
%     r = zeros(nz,1);      %Vector to store the row indexes
    c = zeros(m+1,1);     %Vector to store where column begins
    c(1) = 1;
    
    [r,jj,v] = find(A);
    nz = length(r);
        
    for j=1:nz
        k = jj(j)+1;
        c(k) = c(k)+1;
    end
    for j=2:m+1
        c(j) = c(j)+c(j-1);
    end

end