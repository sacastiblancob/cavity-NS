function [U,V] = f_vel(X,Y,t)
% Function for velocity in X direction (U) and velocity in Y direction (V)
%
% Based on: DG-FTLE: Lagrangian coherent structures with high-order
%           discontinuous-Galerkin methods
%               Daniel A. Nelson, Gustaaf B. Jacobs, 2015

A=0.1;
% A=1;
eps=0.1;
ome=(2*pi)/10;

U = -pi.*A.*sin(pi*(eps*sin(ome*t)*X.*X + (1-2*eps*sin(ome*t)).*X)).*cos(pi.*Y);
V = pi.*A.*cos(pi*(eps*sin(ome*t)*X.*X + (1-2*eps*sin(ome*t)).*X)).*sin(pi.*Y).*...
    (2*eps*sin(ome*t)*X + (1-2*eps*sin(ome*t)));

end