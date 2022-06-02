function [FTLEf,FTLEb,XDGgt2,YDGgt2,XDGgt1,YDGgt1,XDGgt0,YDGgt0] = ...
    FTLE_dgelem_fvel_old(ep,et,XDGgt0,YDGgt0,dt,NT,t1,Nx,Ny,Jx,Jy,XDD,YDD)
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
for t=t1-NT*dt:dt:t1-dt
    [XDGgt1,YDGgt1] = AB3_fvel(XDGgt1,YDGgt1,dt,t);
end

%AB3 for particle tracing from t1 to t2
XDGgt2 = XDGgt1;
YDGgt2 = YDGgt1;
for t=t1:dt:t1+dt*(NT-1)
    [XDGgt2,YDGgt2] = AB3_fvel(XDGgt2,YDGgt2,dt,t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP C - Computing Forward-time FTLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Computing parametrizations
% xb1 = XDGgt1(1,:);   yb1 = YDGgt1(1,:);
% xb2 = XDGgt1(:,end); yb2 = YDGgt1(:,end);
% xb3 = XDGgt1(end,:); yb3 = YDGgt1(end,:);
% xb4 = XDGgt1(:,1);   yb4 = YDGgt1(:,1);
% [ppgx,ppgy,ppgpx,ppgpy] = param(ep,et,xb1,xb2,xb3,xb4,yb1,yb2,yb3,yb4);
% 
% %Mapping to local element
% x1 = [XDGgt1(1,1); YDGgt1(1,1)];
% x2 = [XDGgt1(1,end); YDGgt1(1,end)];
% x3 = [XDGgt1(end,end); YDGgt1(end,end)];
% x4 = [XDGgt1(end,1); YDGgt1(end,1)];
% [XDGlt1,YDGlt1] = map_to_local(ep,et,XDGgt1,YDGgt1,ppgx,ppgy,ppgpx,ppgpy,x1,x2,x3,x4);

%Mapping to local element
[XDGlt1,YDGlt1] = map_to_local_dg(ep,et,XDGgt1,YDGgt1,XDD,YDD);

%Evaluating the interpolation Ipp.^-1 (here PSI = [XDGgt2,YDGgt2] and THETA = [XDGgt0,YDGgt0] is computed
% [XDGgt2,YDGgt2] = interpol(ep,et,XDGlt1,YDGlt1,XDGgt2,YDGgt2);
% [XDGgt0,YDGgt0] = interpol(ep,et,XDGlt1,YDGlt1,XDGgt0,YDGgt0);
XDGgt2 = lagrange_2D_interpolation(ep,et,XDGlt1,YDGlt1,XDGgt2,Nx,Ny);
YDGgt2 = lagrange_2D_interpolation(ep,et,XDGlt1,YDGlt1,YDGgt2,Nx,Ny);
XDGgt0 = lagrange_2D_interpolation(ep,et,XDGlt1,YDGlt1,XDGgt0,Nx,Ny);
YDGgt0 = lagrange_2D_interpolation(ep,et,XDGlt1,YDGlt1,YDGgt0,Nx,Ny);

%Computing Green-Cauchy tensor and FTLE
FTLEf = zeros(size(XDGgt0));
FTLEb = zeros(size(XDGgt0));

% [dxdep,dydep,dxdet,dydet] = derv_epet(Nx,Ny,ep,et,ppgx,ppgy,ppgpx,ppgpy,x1,x2,x3,x4);
[dxdep,dydep,dxdet,dydet] = derv_epet_dg(ep,et,XDGgt1,YDGgt1,XDD,YDD);

dXDGgt2dep = (XDD*XDGgt2')';
dYDGgt2dep = (XDD*YDGgt2')';
dXDGgt2det = YDD*XDGgt2;
dYDGgt2det = YDD*YDGgt2;

dphixdx = (1/Jx)*(dXDGgt2dep.*dydet - dXDGgt2det.*dydep);
dphixdy = (1/Jy)*(dXDGgt2det.*dxdep - dXDGgt2dep.*dxdet);
dphiydx = (1/Jx)*(dYDGgt2dep.*dydet - dYDGgt2det.*dydep);
dphiydy = (1/Jy)*(dYDGgt2det.*dxdep - dYDGgt2dep.*dxdet);

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

dXDGgt2dep = (XDD*XDGgt0')';
dYDGgt2dep = (XDD*YDGgt0')';
dXDGgt2det = YDD*XDGgt0;
dYDGgt2det = YDD*YDGgt0;

dphixdx = (1/Jx)*(dXDGgt2dep.*dydet - dXDGgt2det.*dydep);
dphixdy = (1/Jy)*(dXDGgt2det.*dxdep - dXDGgt2dep.*dxdet);
dphiydx = (1/Jx)*(dYDGgt2dep.*dydet - dYDGgt2det.*dydep);
dphiydy = (1/Jy)*(dYDGgt2det.*dxdep - dYDGgt2dep.*dxdet);

for j=1:Nx+1
    for i=1:Ny+1
        gcT = [dphixdx(i,j) dphiydx(i,j);dphixdy(i,j) dphiydy(i,j)];
        gcT = gcT'*gcT;
        FTLEb(i,j) = (1/(dt*NT))*log(sqrt(eigs(gcT,1)));
    end
end

end