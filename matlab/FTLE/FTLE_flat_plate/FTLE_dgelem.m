function [FTLEf,FTLEb,XDGgt2,YDGgt2,XDGgt1,YDGgt1,XDGgt0,YDGgt0] = ...
    FTLE_dgelem(X,Y,U,V,ep,et,XDGgt0,YDGgt0,dt,NT,Nx,Ny,XDD,YDD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes de Finite Time Lyapunov Exponents fields for a
% Discontinuous Galerkin element.
% Velocity field, for instance, comes from a Finite Difference model, or a
% structured mesh.
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
% CFL = 2; NT=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP A - Particle tracing algorithm - Adam-Bashforth 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%AB3 for particle tracing from t0 to t1
XDGgt1 = XDGgt0;
YDGgt1 = YDGgt0;
for t=1:NT
    [XDGgt1,YDGgt1] = AB3(X,Y,U,V,XDGgt1,YDGgt1,dt);
end

%AB3 for particle tracing from t1 to t2
XDGgt2 = XDGgt1;
YDGgt2 = YDGgt1;
for t=1:NT
    [XDGgt2,YDGgt2] = AB3(X,Y,U,V,XDGgt2,YDGgt2,dt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP C - Computing Forward-time FTLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Mapping to local element
[XDGlt1,YDGlt1] = map_to_local_dg(ep,et,XDGgt1,YDGgt1,XDD,YDD);

%Evaluating the interpolation Ipp.^-1 (here PSI = [XDGgt2,YDGgt2] and THETA = [XDGgt0,YDGgt0] is computed
XDGgt2 = lagrange_2D_interpolation(ep,et,XDGlt1,YDGlt1,XDGgt2,Nx,Ny);
YDGgt2 = lagrange_2D_interpolation(ep,et,XDGlt1,YDGlt1,YDGgt2,Nx,Ny);
XDGgt0 = lagrange_2D_interpolation(ep,et,XDGlt1,YDGlt1,XDGgt0,Nx,Ny);
YDGgt0 = lagrange_2D_interpolation(ep,et,XDGlt1,YDGlt1,YDGgt0,Nx,Ny);

%Computing Green-Cauchy tensor and FTLE
FTLEf = zeros(size(XDGgt0));
FTLEb = zeros(size(XDGgt0));

[dxdep,dydep,dxdet,dydet] = derv_epet_dg(ep,et,XDGgt1,YDGgt1,XDD,YDD);

%Computing Mapping Jacobian
J = dxdep.*dydet - dydep.*dxdet;

%Computing derivatives of PHI with respect to local element
dXDGgt2dep = (XDD*XDGgt2')';
dYDGgt2dep = (XDD*YDGgt2')';
dXDGgt2det = YDD*XDGgt2;
dYDGgt2det = YDD*YDGgt2;

%Computing derivatives with respect to global element
dphixdx = (1./J).*(dXDGgt2dep.*dydet - dXDGgt2det.*dydep);
dphixdy = (1./J).*(dXDGgt2det.*dxdep - dXDGgt2dep.*dxdet);
dphiydx = (1./J).*(dYDGgt2dep.*dydet - dYDGgt2det.*dydep);
dphiydy = (1./J).*(dYDGgt2det.*dxdep - dYDGgt2dep.*dxdet);

for j=1:Nx+1
    for i=1:Ny+1
        gcT = [dphixdx(i,j) dphiydx(i,j);dphixdy(i,j) dphiydy(i,j)];
        gcT = gcT'*gcT;
        FTLEf(i,j) = (1/(dt*NT))*log(sqrt(eigs(gcT,1)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP D - Computing Backward-time FTLE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Computing derivatives of OME with respect to local element
dXDGgt2dep = (XDD*XDGgt0')';
dYDGgt2dep = (XDD*YDGgt0')';
dXDGgt2det = YDD*XDGgt0;
dYDGgt2det = YDD*YDGgt0;

%Computing derivatives with respect to global element
dphixdx = (1./J).*(dXDGgt2dep.*dydet - dXDGgt2det.*dydep);
dphixdy = (1./J).*(dXDGgt2det.*dxdep - dXDGgt2dep.*dxdet);
dphiydx = (1./J).*(dYDGgt2dep.*dydet - dYDGgt2det.*dydep);
dphiydy = (1./J).*(dYDGgt2det.*dxdep - dYDGgt2dep.*dxdet);

for j=1:Nx+1
    for i=1:Ny+1
        gcT = [dphixdx(i,j) dphiydx(i,j);dphixdy(i,j) dphiydy(i,j)];
        gcT = gcT'*gcT;
        FTLEb(i,j) = (1/(dt*NT))*log(sqrt(eigs(gcT,1)));
    end
end


end