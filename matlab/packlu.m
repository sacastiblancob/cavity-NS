function [LU] = packlu(A)
%Computes de LU decompostion for matrix A and stores it in a packed way,
%thats it, into the same matrix LU, where the diagonal and the triangular
%superior part correspond with the matrix U, and the inferior part
%correspond with the entries of matrix L

n = length(A);

LU = A;

for j = 1:n-1
    % Compute a column of L, then update submatrix
    LU(j+1:n,j) = LU(j+1:n,j)/LU(j,j);
    LU(j+1:n,j+1:n) = LU(j+1:n,j+1:n)-LU(j+1:n,j)*LU(j,j+1:n);
end

end