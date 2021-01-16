function x = solpacklu(LU,b)
%Computes the forward substitution and backward substitution for compute
%the solution to the system LUx=b, with LU being a packed-LU decomposition
%matrix

n = length(LU);

x = b;

for i = 1:n
    x(i) = x(i)-LU(i,1:i-1)*x(1:i-1);
end
for i = n:-1:1
    x(i) = ( x(i)-LU(i,i+1:n)*x(i+1:n) )/LU(i,i);
end

end