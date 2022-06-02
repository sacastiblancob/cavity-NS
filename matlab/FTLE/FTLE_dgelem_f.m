function [FTLEf,XDGgt1,YDGgt1] = ...
    FTLE_dgelem_f(X,Y,U,V,XDGgt0,YDGgt0,dt,NT,Nx,Ny,Jx,Jy,XDD,YDD)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP C - Computing Forward-time FTLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Computing Green-Cauchy tensor and FTLE
FTLEf = zeros(size(XDGgt0));

dXDGgt1dep = (XDD*XDGgt1')';
dYDGgt1dep = (XDD*YDGgt1')';
dXDGgt1det = YDD*XDGgt1;
dYDGgt1det = YDD*YDGgt1;

dphixdx = (1/Jx)*(dXDGgt1dep - dXDGgt1det);
dphixdy = (1/Jy)*(dXDGgt1det - dXDGgt1dep);
dphiydx = (1/Jx)*(dYDGgt1dep - dYDGgt1det);
dphiydy = (1/Jy)*(dYDGgt1det - dYDGgt1dep);

for j=1:Nx+1
    for i=1:Ny+1
        gcT = [dphixdx(i,j) dphiydx(i,j);dphixdy(i,j) dphiydy(i,j)];
        gcT = gcT'*gcT;
        FTLEf(i,j) = (1/(dt*NT))*log(sqrt(eigs(gcT,1)));
    end
end

end