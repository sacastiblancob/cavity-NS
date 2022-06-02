function [Kt0,K0,dxdep,dydep,XDGgt1,YDGgt1] = ...
    lagcurvchange_dgelem(X,Y,U,V,ep,et,XDGgt0,YDGgt0,dt,NT,NT1,Nx,Ny,XDD,YDD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes de Lagrangian Curvature Change (kt0) with DG operators
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

%%%%%%%%%%%%%%5
% FIRST, computing initial curvature scalars

%AB3 for particle tracing from parallel material lines to t0
XDGgt1 = XDGgt0;
YDGgt1 = YDGgt0;
for t=1:NT
    [XDGgt1,YDGgt1] = AB3(X,Y,U,V,XDGgt1,YDGgt1,dt);
end

%Computing Curvature Scalars
K0 = zeros(size(XDGgt0));

[dxdep,dydep] = derv_epet_dg(ep,et,XDGgt1,YDGgt1,XDD,YDD);
% dxdep = -dxdep; dydep = -dydep;
% [dxdep2,dydep2] = derv_epet_dg(ep,et,dxdep,dydep,XDD,YDD);
dxdep2 = (XDD*dxdep')';
dydep2 = (XDD*dydep')';
dxdep2 = -dxdep2; dydep2 = -dydep2;
% dxdep2 = YDD*dxdep;
% dydep2 = YDD*dydep;

for j=1:Nx+1
    for i=1:Ny+1
        K0(i,j) = (dxdep2(i,j)*dydep(i,j) - dydep2(i,j)*dxdep(i,j))/...
            (dxdep(i,j)^2 + dydep(i,j)^2)^(1/3);
    end
end

%%%%%%%%%%%%%%5
% SECOND, computing Lagrangian Curvature Change for t1 from NT to NT+NT1

I1 = zeros(size(XDGgt0));
I2 = zeros(size(XDGgt0));

for t=1:NT1

    %Computing tangent vectors
    [dxdep,dydep,dxdet,dydet] = derv_epet_dg(ep,et,XDGgt1,YDGgt1,XDD,YDD);
    % dxdep = -dxdep; dydep = -dydep;
    
    %Computing Mapping Jacobian
    J = dxdep.*dydet - dydep.*dxdet;
    
    %Interpolating velocity field
    UDG1 = interp2(X,Y,U,XDGgt1,YDGgt1,'makima');
    VDG1 = interp2(X,Y,V,XDGgt1,YDGgt1,'makima');
    
    %Computing du/dep du/det dv/dep dv/det
    dudep = (XDD*UDG1')'; dudet = YDD*UDG1;
    dvdep = (XDD*VDG1')'; dvdet = YDD*VDG1;
    
    %Computing du/dx du/dy dv/dx dv/dy
    dudx = (1./J).*(dudep.*dydet - dudet.*dydep);
    dudy = (1./J).*(dudet.*dxdep - dudep.*dxdet);
    dvdx = (1./J).*(dvdep.*dydet - dvdet.*dydep);
    dvdy = (1./J).*(dvdet.*dxdep - dvdep.*dxdet);
    
    %Computing d2u/dx2 d2u/dy2 d2u/dxdy d2v/dxdy d2v/dx2 d2v/dy2
    dudxdep = (XDD*dudx')'; dudxdet = YDD*dudx;
    dvdxdep = (XDD*dvdx')'; dvdxdet = YDD*dvdx;
    dudydep = (XDD*dudy')'; dudydet = YDD*dudy;
    dvdydep = (XDD*dvdy')'; dvdydet = YDD*dvdy;
    
    d2udx2 = (1./J).*(dudxdep.*dydet - dudxdet.*dydep);
    d2udxdy = (1./J).*(dudxdet.*dxdep - dudxdep.*dxdet);
    d2udy2 = (1./J).*(dudydet.*dxdep - dudydep.*dxdet);
    d2vdxdy = (1./J).*(dvdydep.*dydet - dvdydet.*dydep);
    d2vdy2 = (1./J).*(dvdydet.*dxdep - dvdydep.*dxdet);
    d2vdx2 = (1./J).*(dvdxdep.*dydet - dvdxdet.*dydep);
    
    %computing nabla*omega (gradient of vorticity)
    dodep = (XDD*(dvdx - dudy)')'; dodet = YDD*(dvdx - dudy);
    
    dodx = (1./J).*(dodep.*dydet - dodet.*dydep);
    dody = (1./J).*(dodet.*dxdep - dodep.*dxdet);
    
    %Computing terms for integrating kt0
    nrr = dxdep.^2 + dxdep.^2;
    rs = dxdep.*(dudx.*dxdep + 0.5*(dudy+dvdx).*dydep) + dydep.*(0.5*(dudy+dvdx).*dxdep + dvdy.*dydep);
    rs = rs./nrr;
    rns = dydep.*((d2udx2.*dxdep + d2udxdy.*dydep).*dxdep + (0.5*(d2udxdy + d2vdx2).*dxdep + 0.5*(d2udy2 + d2vdxdy).*dydep).*dydep) -...
        dxdep.*((0.5*(d2udxdy+d2vdx2).*dxdep + 0.5*(d2udy2+d2vdxdy).*dydep).*dxdep + (d2vdxdy.*dxdep + d2vdy2.*dydep).*dydep);
    rns = rns./(nrr.^(3/2));
    nor = dodx.*dxdep + dody.*dydep;
    nor = nor./(2*(nrr.^(1/2)));

    %Computing integrals with euler
    I1 = I1 + dt*rs;
    I2 = I2 + dt*((rns-nor).*exp(3*I1));

    %Advecting particles
    [XDGgt1,YDGgt1] = AB3(X,Y,U,V,XDGgt1,YDGgt1,dt);

end

Kt0 = exp(-3.*I1).*I2 - K0;

end
