function [lk] = lagrange(xk,xi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the lagrange polynomial based on nodes xi, on
% values estored in xk for every i
% i.e., computes Vandermonde matrix based on lagrange polynomials
%  if xk==xi, lk = I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(xk);
lk = zeros(n,n);
for k=1:n
    for i=1:n
        num = xk(k) - xi;
        den = xi(i) - xi;
        lk(i,k) = prod(num(den~=0)./den(den~=0));
    end
end

end


