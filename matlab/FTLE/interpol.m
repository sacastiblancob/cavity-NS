function [XDGg2,YDGg2] = interpol(ep,et,XDGl1,YDGl1,XDGg1,YDGg1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the interpolation from advected nodes

%
% Based on: DG-FTLE: Lagrangian coherent structures with high-order
%           discontinuous-Galerkin methods
%               Daniel A. Nelson, Gustaaf B. Jacobs, 2015

n = length(ep);
XDGg2 = zeros(n,n);
YDGg2 = zeros(n,n);

for i=1:n
%     li = lagrange(ep,XDGl1(i,:));
%     mj = lagrange(et,YDGl1(:,i));
%     Ii = li.*mj;
    li = lagrange(ep,XDGl1(i,:));
    mj = lagrange(et,YDGl1(:,i));
    Ii = (mj*li)';
    XDGg2(:,i) = Ii*XDGg1(:,i);
    YDGg2(:,i) = Ii*YDGg1(:,i);

%     XDGg2(:,i) = li*(mj*XDGg1(:,i));
%     YDGg2(:,i) = li*(mj*YDGg1(:,i));
end

end

