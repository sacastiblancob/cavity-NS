function [XDGlg1,YDGlg1] = map_to_global_dg(ep,et,XDGg1,YDGg1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function map local coordinates to global coordinates
%YDGgt1
%                       b3(xb3,yb3)           
%            ---------------------------------
%            |                               |
%            |                               |
% b4(xb4,yb4)|                               |b2(xb2,yb2)
%            |                               |
%            |                               |
%            ---------------------------------
%                       b1(xb1,yb1)
%
%
% Entries:
%  - ep, et: local-master element coordinates
%  - XDGg1, YDGg1 : global coordinates points for mapping into local
%  element
%  - ppgx, ppgy : parametrizations for boundaries of the element (see param.m)
%  - ppgpx, ppgpy : derivatives of parametrizations for boundaries of the element (see param.m)
%  - x1,x2,x3,x4 : corners coordinates of deformed element, column vector
%                   x1 = [x;y]_1 ; x2 = [x;y]_2 ; etc...
%
% Based on: DG-FTLE: Lagrangian coherent structures with high-order
%           discontinuous-Galerkin methods
%               Daniel A. Nelson, Gustaaf B. Jacobs, 2015
%
% Sergio Castiblanco-Ballesteros
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Extracting parametrizations of the boundaries of element
xb1 = XDGg1(1,:)';   yb1 = YDGg1(1,:)';
xb2 = XDGg1(:,end); yb2 = YDGg1(:,end);
xb3 = XDGg1(end,:)'; yb3 = YDGg1(end,:)';
xb4 = XDGg1(:,1);   yb4 = YDGg1(:,1);

%Corners of the deformed element
x1 = [XDGg1(1,1); YDGg1(1,1)];
x2 = [XDGg1(1,end); YDGg1(1,end)];
x3 = [XDGg1(end,end); YDGg1(end,end)];
x4 = [XDGg1(end,1); YDGg1(end,1)];

%Size of the element
sizex = size(XDGg1);
Nx = sizex(2)-1;      %Polynomial order in x-direction
Ny = sizex(1)-1;      %Polynomial order in y-direction

%Results will be stored here
XDGlg1 = XDGg1;
YDGlg1 = YDGg1;

%GOING INTO NEWTON'S METHOD
for j=2:Nx
    for i=2:Ny

        %Evaluating parametrizations
        g1 = [xb1(j); yb1(j)];
        g2 = [xb2(i); yb2(i)];
        g3 = [xb3(j); yb3(j)];
        g4 = [xb4(i); yb4(i)];

        ep0 = ep(j);
        et0 = et(i);

        %x(ep,et) evaluated at first guess
        F1 = (1-et0)*g1 + et0*g3 + (1-ep0)*g4 + ep0*g2 - x1*(1-ep0)*(1-et0) - ...
            x2*ep0*(1-et0) - x3*ep0*et0 - x4*(1-ep0)*et0;

        XDGlg1(i,j) = F1(1);
        YDGlg1(i,j) = F1(2);
    end
end


end