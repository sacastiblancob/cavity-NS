function [dxdep,dydep,dxdet,dydet] = derv_epet(Nx,Ny,ep,et,ppgx,ppgy,ppgpx,ppgpy,x1,x2,x3,x4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes de Derivatives with respect to local-master
% coordinates epsilon (ep) eta (et).
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

dxdep = zeros(Ny+1,Nx+1);
dydep = zeros(Ny+1,Nx+1);
dxdet = zeros(Ny+1,Nx+1);
dydet = zeros(Ny+1,Nx+1);

for j=1:Nx+1
    for i=1:Ny+1

        ep0 = ep(j);
        et0 = et(i);

        %Evaluating parametrizations
        g1 = [ppval(ppgx(1),ep0); ppval(ppgy(1),ep0)];
        g2 = [ppval(ppgx(2),et0); ppval(ppgy(2),et0)];
        g3 = [ppval(ppgx(3),ep0); ppval(ppgy(3),ep0)];
        g4 = [ppval(ppgx(4),et0); ppval(ppgy(4),et0)];

        %Evuating derivations of parametrizations
        g1p = [ppval(ppgpx(1),ep0); ppval(ppgpy(1),ep0)];
        g2p = [ppval(ppgpx(2),et0); ppval(ppgpy(2),et0)];
        g3p = [ppval(ppgpx(3),ep0); ppval(ppgpy(3),ep0)];
        g4p = [ppval(ppgpx(4),et0); ppval(ppgpy(4),et0)];

        %Computing the Jacobian
        dxydep = x1*(1-et0) - x2*(1-et0) - x3*et0 + x4*et0 + g1p*(1-et0) + g3p*et0 + g2 - g4;
        dxydet = x1*(1-ep0) + x2*ep0 - x3*ep0 - x4*(1-ep0) + g2p*ep0 + g4p*(1-ep0) - g1 + g3;

        dxdep(i,j) = dxydep(1);
        dydep(i,j) = dxydep(2);
        dxdet(i,j) = dxydet(1);
        dydet(i,j) = dxydet(2);
    end
end

end