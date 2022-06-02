function [le] = lagrange_point(xe,xi,i)
% This function evaluates lagrange function i, with elements in xi as basis
% in the point xe

num = xe - xi;
den = xi(i) - xi;
le = prod(num(1:length(xi)~=i)./den(1:length(xi)~=i));

end