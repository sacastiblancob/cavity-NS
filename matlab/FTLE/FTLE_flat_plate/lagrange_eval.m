function [le] = lagrange_eval(xe,xi,yi)
% This function evaluates lagrange polynomial with elements in xi as basis
% in the point xe

le = 0;
for i=1:length(xi)
    num = xe - xi;
    den = xi(i) - xi;
    le = le + yi(i)*prod(num(den~=0)./den(den~=0));
end

end