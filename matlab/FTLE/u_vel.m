function [U,V] = f_vel(X,Y,t,A,eps,ome)
% Function for velocity in X direction
%
% Based on: DG-FTLE: Lagrangian coherent structures with high-order
%           discontinuous-Galerkin methods
%               Daniel A. Nelson, Gustaaf B. Jacobs, 2015

U = -pi*A*sin(pi*(eps*sin(ome*t)*X.*X + (1-2*eps*sin(ome*t)).*X))*cos(pi.*Y);

end