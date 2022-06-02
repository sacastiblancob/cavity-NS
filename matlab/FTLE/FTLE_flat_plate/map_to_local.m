function [XDGl1,YDGl1] = map_to_local(ep,et,XDGg1,YDGg1,ppgx,ppgy,ppgpx,ppgpy,x1,x2,x3,x4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function map global coordinates to local coordinates using
% Newton-Rhapson method
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sizex = size(XDGg1);
Nx = sizex(2)-1;      %Polynomial order in x-direction
Ny = sizex(1)-1;      %Polynomial order in y-direction

%Results will be stored here
XDGl1 = kron(ones(Ny+1,1),ep);
YDGl1 = kron(et,ones(1,Nx+1));

%GOING INTO NEWTON'S METHOD
for j=2:Nx
    for i=2:Ny

        %Global coordinates of the point
        xy0 = [XDGg1(i,j); YDGg1(i,j)];

        %First guess
        epet0 = [ep(j) et(i)];

        for k=1:100
            %initial guess
            ep0 = epet0(1);
            et0 = epet0(2);

            %Evaluating parametrizations
            g1 = [ppval(ppgx(1),ep0); ppval(ppgy(1),ep0)];
            g2 = [ppval(ppgx(2),et0); ppval(ppgy(2),et0)];
            g3 = [ppval(ppgx(3),ep0); ppval(ppgy(3),ep0)];
            g4 = [ppval(ppgx(4),et0); ppval(ppgy(4),et0)];

            %Evuating derivatives of parametrizations
            g1p = [ppval(ppgpx(1),ep0); ppval(ppgpy(1),ep0)];
            g2p = [ppval(ppgpx(2),et0); ppval(ppgpy(2),et0)];
            g3p = [ppval(ppgpx(3),ep0); ppval(ppgpy(3),ep0)];
            g4p = [ppval(ppgpx(4),et0); ppval(ppgpy(4),et0)];

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
            if norm(abs(epet1 - epet0)) < 1E-8
                break
            end

            %updating
            epet0 = epet1;
        end
%         k
        XDGl1(i,j) = epet1(1);
        YDGl1(i,j) = epet1(2);
    end
end

end