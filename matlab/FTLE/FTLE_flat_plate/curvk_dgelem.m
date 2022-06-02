function [K0,dxdep,dydep,XDGgt1,YDGgt1,dxdep2,dydep2] = ...
    curvk_dgelem(X,Y,U,V,ep,et,XDGgt0,YDGgt0,dt,NT,Nx,Ny,XDD,YDD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes de Curvature Scalars (k) with DG operators
%
% Inputs:
%   - X,Y : matrices with X and Y coordinates of the velocity field.
%   - U,V : matrices with x-component and y-component of velocity field.
%   - ep,et : local-master unitary element coordinates vectors (ep(epsilon-x), et(eta-y))
%   - XDGg, YDGg : global element coordinates matrices.
%   - CFL : Courant Number for defining convective-time units
%   - NT : number of convective-time units for particle tracing

%
% Based on: DG-FTLE: Lagrangian coherent structures with high-order
%           discontinuous-Galerkin methods
%               Daniel A. Nelson, Gustaaf B. Jacobs, 2015
%         : Exact theory of material spike formation in flow separation.
%               Serra, Vétel, Haller, 2018
% CFL = 2; NT=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%AB3 for particle tracing from parallel material lines to t0
XDGgt1 = XDGgt0;
YDGgt1 = YDGgt0;
for t=1:NT
    [XDGgt1,YDGgt1] = AB3(X,Y,U,V,XDGgt1,YDGgt1,dt);
end

% %Computing Curvature Scalars
% K0 = zeros(size(XDGgt0));
% 
% [dxdep,dydep] = derv_epet_dg(ep,et,XDGgt1,YDGgt1,XDD,YDD);
% dxdep = -dxdep; dydep = -dydep;
% % [dxdep2,dydep2] = derv_epet_dg(ep,et,dxdep,dydep,XDD,YDD);
% dxdep2 = (XDD*dxdep')';
% dydep2 = (XDD*dydep')';
% % dxdep2 = (1/Jx)*(XDD*dxdep')';
% % dydep2 = (1/Jx)*(XDD*dydep')';
% dxdep2 = -dxdep2; dydep2 = -dydep2;
% % dxdep2 = YDD*dxdep;
% % dydep2 = YDD*dydep;
% 
% for j=1:Nx+1
%     for i=1:Ny+1
%         K0(i,j) = (dxdep2(i,j)*dydep(i,j) - dydep2(i,j)*dxdep(i,j))/...
%             (dxdep(i,j)^2 + dydep(i,j)^2)^(1/3);
%     end
% end

%Computing Curvature Scalars
K0 = zeros(size(XDGgt0));

[dxdep,dydep,dxdet,dydet] = derv_epet_dg(ep,et,XDGgt1,YDGgt1,XDD,YDD);

%Computing Mapping Jacobian
J = dxdep.*dydet - dydep.*dxdet;

%Computing d2x/dep2 d2y/dep2
dxdep = -dxdep; dydep = -dydep;
d2xdep2 = (XDD*dxdep')'; d2xdepdet = YDD*dxdep;
d2ydep2 = (XDD*dydep')'; d2ydepdet = YDD*dydep;

dxdep2 = (1./J).*(d2xdep2.*dydet - d2xdepdet.*dydep);
dydep2 = (1./J).*(d2ydep2.*dydet - d2ydepdet.*dydep);
dxdep2 = -dxdep2; dydep2 = -dydep2;

for j=1:Nx+1
    for i=1:Ny+1
        K0(i,j) = (dxdep2(i,j)*dydep(i,j) - dydep2(i,j)*dxdep(i,j))/...
            (dxdep(i,j)^2 + dydep(i,j)^2)^(1/3);
    end
end

end
