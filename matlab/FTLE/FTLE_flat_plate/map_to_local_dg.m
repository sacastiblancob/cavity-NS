function [XDGl1,YDGl1] = map_to_local_dg(ep,et,XDGg1,YDGg1,XDD,YDD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function map global coordinates to local coordinates using
% Newton-Rhapson method
%
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

%Results will be stored here
XDGl1 = kron(ones(Ny+1,1),ep);
YDGl1 = kron(et,ones(1,Nx+1));

%GOING INTO NEWTON'S METHOD
for j=2:Nx
    for i=2:Ny

        %Evaluating parametrizations
        g1 = [xb1(j); yb1(j)];
        g2 = [xb2(i); yb2(i)];
        g3 = [xb3(j); yb3(j)];
        g4 = [xb4(i); yb4(i)];

        %Evaluating derivations of parametrizations
        g1p = [xb1p(j); yb1p(j)];
        g2p = [xb2p(i); yb2p(i)];
        g3p = [xb3p(j); yb3p(j)];
        g4p = [xb4p(i); yb4p(i)];

        %Global coordinates of the point
        xy0 = [XDGg1(i,j); YDGg1(i,j)];

        %First guess vector
        epet0 = [ep(j);et(i)];

        %initial guess
        ep0 = epet0(1);
        et0 = epet0(2);

        for k=1:10 %Up to 10 iterations, slow convergence!
            
            %Computing the Jacobian
            dxdep = x1*(1-et0) - x2*(1-et0) - x3*et0 + x4*et0 + g1p*(1-et0) + g3p*et0 + g2 - g4;
            dxdet = x1*(1-ep0) + x2*ep0 - x3*ep0 - x4*(1-ep0) + g2p*ep0 + g4p*(1-ep0) - g1 + g3;
            J = [dxdep dxdet];

            %x(ep,et) evaluated at first guess
            F1 = (1-et0)*g1 + et0*g3 + (1-ep0)*g4 + ep0*g2 - x1*(1-ep0)*(1-et0) - ...
                x2*ep0*(1-et0) - x3*ep0*et0 - x4*(1-ep0)*et0;

            %Next guess
            epet1 = epet0 - J\(F1-xy0);

            %Comparing norms
%             k
%             abs(epet1 - epet0)  
            if abs(epet1 - epet0) < 5E-3
                break
            end

            %Updating
            epet0 = epet1;
            ep0 = epet0(1);
            et0 = epet0(2);

            %Evaluating parametrizations
            g1 = [lagrange_eval(ep0,ep,xb1); lagrange_eval(ep0,ep,yb1)];
            g2 = [lagrange_eval(et0,et,xb2); lagrange_eval(et0,et,yb2)];
            g3 = [lagrange_eval(ep0,ep,xb3); lagrange_eval(ep0,ep,yb3)];
            g4 = [lagrange_eval(et0,et,xb4); lagrange_eval(et0,et,yb4)];

            %Evuating derivatives of parametrizations
            g1p = [lagrange_eval(ep0,ep,xb1p); lagrange_eval(ep0,ep,yb1p)];
            g2p = [lagrange_eval(et0,et,xb2p); lagrange_eval(et0,et,yb2p)];
            g3p = [lagrange_eval(ep0,ep,xb3p); lagrange_eval(ep0,ep,yb3p)];
            g4p = [lagrange_eval(et0,et,xb4p); lagrange_eval(et0,et,yb4p)];
        end
%         k
        XDGl1(i,j) = epet1(1);
        YDGl1(i,j) = epet1(2);
    end
end

end