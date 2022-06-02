function [dxdep,dydep,dxdet,dydet] = derv_epet_dg(ep,et,XDGg1,YDGg1,XDD,YDD)
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

%Extracting parametrizations of the boundaries of element
xb1 = XDGg1(1,:)';   yb1 = YDGg1(1,:)';
xb2 = XDGg1(:,end); yb2 = YDGg1(:,end);
xb3 = XDGg1(end,:)'; yb3 = YDGg1(end,:)';
xb4 = XDGg1(:,1);   yb4 = YDGg1(:,1);

%Computing derivatives of the boundaries of element
xb1p = XDD*xb1; yb1p = XDD*yb1;
xb2p = YDD*xb2; yb2p = YDD*yb2;
xb3p = XDD*xb3; yb3p = XDD*yb3;
xb4p = YDD*xb4; yb4p = YDD*yb4;

%Corners of the deformed element
x1 = [XDGg1(1,1); YDGg1(1,1)];
x2 = [XDGg1(1,end); YDGg1(1,end)];
x3 = [XDGg1(end,end); YDGg1(end,end)];
x4 = [XDGg1(end,1); YDGg1(end,1)];

%Size of the element
sizex = size(XDGg1);
Nx = sizex(2)-1;      %Polynomial order in x-direction
Ny = sizex(1)-1;      %Polynomial order in y-direction

% %Results will be stored here
% dxdep = zeros(Ny+1,Nx+1);
% dydep = zeros(Ny+1,Nx+1);
% dxdet = zeros(Ny+1,Nx+1);
% dydet = zeros(Ny+1,Nx+1);
% 
% for j=1:Nx+1
%     for i=1:Ny+1
% 
%         ep0 = ep(j);
%         et0 = et(i);
% 
%         %Evaluating parametrizations
%         g1 = [xb1(j); yb1(j)];
%         g2 = [xb2(i); yb2(i)];
%         g3 = [xb3(j); yb3(j)];
%         g4 = [xb4(i); yb4(i)];
% 
%         %Evaluating derivations of parametrizations
%         g1p = [xb1p(j); yb1p(j)];
%         g2p = [xb2p(i); yb2p(i)];
%         g3p = [xb3p(j); yb3p(j)];
%         g4p = [xb4p(i); yb4p(i)];
% 
%         %Computing the Jacobian
%         dxydep = x1*(1-et0) - x2*(1-et0) - x3*et0 + x4*et0 + g1p*(1-et0) + g3p*et0 + g2 - g4;
%         dxydet = x1*(1-ep0) + x2*ep0 - x3*ep0 - x4*(1-ep0) + g2p*ep0 + g4p*(1-ep0) - g1 + g3;
% 
%         dxdep(i,j) = dxydep(1);
%         dydep(i,j) = dxydep(2);
%         dxdet(i,j) = dxydet(1);
%         dydet(i,j) = dxydet(2);
%     end
% end

% dxdepv = x1(1)*(1-et) - x2(1)*(1-et) - x3(1)*et + x4(1)*et + xb2 - xb4;
% dxdepm = (1-et)*xb1p' + et*xb3p';
% dxdep1 = dxdepm + kron(dxdepv,ones(1,Nx+1));
dxdep = (1-et)*xb1p' + et*xb3p' + kron(x1(1)*(1-et) - x2(1)*(1-et) - x3(1)*et + x4(1)*et + xb2 - xb4,ones(1,Nx+1));

% dydepv = x1(2)*(1-et) - x2(2)*(1-et) - x3(2)*et + x4(2)*et + yb2 - yb4;
% dydepm = (1-et)*yb1p' + et*yb3p';
% dydep1 = dydepm + kron(dydepv,ones(1,Nx+1));
dydep = (1-et)*yb1p' + et*yb3p' + kron(x1(2)*(1-et) - x2(2)*(1-et) - x3(2)*et + x4(2)*et + yb2 - yb4,ones(1,Nx+1));

% dxdetv = x1(1)*(1-ep) + x2(1)*ep - x3(1)*ep - x4(1)*(1-ep) - xb1' + xb3';
% dxdetm = xb2p*ep + xb4p*(1-ep);
% dxdet1 = dxdetm + kron(ones(Ny+1,1),dxdetv);
dxdet = xb2p*ep + xb4p*(1-ep) + kron(ones(Ny+1,1),x1(1)*(1-ep) + x2(1)*ep - x3(1)*ep - x4(1)*(1-ep) - xb1' + xb3');

% dydetv = x1(2)*(1-ep) + x2(2)*ep - x3(2)*ep - x4(2)*(1-ep) - yb1' + yb3';
% dydetm = yb2p*ep + yb4p*(1-ep);
% dydet1 = dydetm + kron(ones(Ny+1,1),dydetv);
dydet = yb2p*ep + yb4p*(1-ep) + kron(ones(Ny+1,1),x1(2)*(1-ep) + x2(2)*ep - x3(2)*ep - x4(2)*(1-ep) - yb1' + yb3');



end










