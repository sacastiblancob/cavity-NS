function [le] = lagrange_point2(xe,xi)
% This function evaluates lagrange polynomial with elements in xi as basis
% in the point xe

le = 0;
for i=1:length(xi)
    num = xe - xi;
    den = xi(i) - xi;
    le = le + prod(num(den~=0)./den(den~=0));
end

end