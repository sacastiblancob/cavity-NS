function x = csc_solpacklu(LUv, LUr, LUc, b)
% This function compute the forward and backward sustitution for the system
% LUx = b, where LU should be entered in pack_csc notation, that's it:
% packLU in Compressed Sparse Column notation.
% Enries are:
%       LUv = vector with values of packLU matrix in CSC format
%       LUr = vector with row indexes of packLU matrix in CSC format
%       LUc = vector with column indexes of packLU matrix in CSC format

x = b;
m = length(LUc)-LUc(1);

%forward sustitution
for j=1:m
    for i = LUc(j):LUc(j+1)-LUc(1)
        if j<LUr(i)
            x(LUr(i)) = x(LUr(i)) - LUv(i)*x(j);
        end
    end
end

%backward sustitution
for j=m:-1:1
    for i = LUc(j+1)-LUc(1):-1:LUc(j)
        if j==LUr(i)
            x(j) = x(j)/LUv(i);
        end
        if j>LUr(i)
            %y(LUr(i)) = y(LUr(i)) - LUv(i)*y(j);
            x(LUr(i)) = x(LUr(i)) - LUv(i)*x(j);
            %LUv(i)
        end
    end
end

end