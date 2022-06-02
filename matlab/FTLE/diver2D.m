function [dudx,dudy,dvdx,dvdy] = diver2D(Uo,Vo,dx,dy,n,m,bound,...
        upboundin,doboundin,leboundin,riboundin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the terms du/dx, du/dy, dv/dx, dv/dy with upwind
% second order scheme (if posible)
%
% Nonlineal advection equation
%
%   du/dt + u du/dx + v du/dy = 0
%
%   dv/dt + u dv/dx + v dv/dy = 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entries
%   Uo, vector with velocities in X direction
%   Uv, vector with velocities in Y direction
%   dx, differential in X
%   dy, differential in Y
%   n, number of nodes in X direction
%   m, number of nodes in Y direction
%   bound, indices of the boundary nodes (b)
%   upboundin, indices with the top internal boundary (u)
%   doboundin, indices with the bottom internal boundary (d)
%   leboundin, indices with the left internal boundary (l)
%   riboundin, indices with the right internal boundary (r)
%
% Grid and what bound, upboundin, doboundin, leboundin, and riboundin means
%
% Enumartion used in the grid
%
%    25   26   27   28   29   30
%
%    19   20   21   22   23   24
%
%    13   14   15   16   17   18
%
%    7    8    9    10   11   12
%
%    1    2    3    4    5    6
%
% bound(b), upboundin(u), doboundin(d), leboundin(l), riboundin(r),
% internal nodes(i)
%
%    b    b    b    b    b    b
%
%    b    ul   u    u    ur   b
%
%    b    l    i    i    r    b
%
%    b    dl   d    d    dr   b
%
%    b    b    b    b    b    b
%
% ul --> top-left internal
% dl --> bottom-left internal
% dr --> bottom-right internal
% ur --> top-right internal
%
% The results in the boundary are set to zero, because they are no needed
% to solve non-lineal advection equation (due to the boundary conditions)
% In other words, the divergence is only computed in the internal nodes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %Allocating memory for the result
    dudx = zeros(n*m,1);
    dudy = zeros(n*m,1);
    dvdx = zeros(n*m,1);
    dvdy = zeros(n*m,1);
    
    %
    %Derivatives in the internal nodes
    %
    for i=1:n*m
        if Uo(i) > 0 && Vo(i) > 0  %which means backward upwind in both directions
            %[doboundin, bottom-right internal node], second order in x, first order in y
            %[leboundin, top-left internal node], first order in x, second order in y
            %[bottom left internal node], first order in both directions
            %[riboundin,upboundin,top-right internal node, pure internal nodes], second order in x and y
            if ismember(i,bound)
                continue
            elseif (ismember(i,doboundin) || i==(2*n-1)) %doboundin, bottom-right internal node
                 dudx(i) = (1/(2*dx))*(3*Uo(i) - 4*Uo(i-1) + Uo(i-2));
                 dudy(i) = (1/dy)*(Uo(i) - Uo(i-n));
                 dvdx(i) = (1/(2*dx))*(3*Vo(i) - 4*Vo(i-1) + Vo(i-2));
                 dvdy(i) = (1/dy)*(Vo(i) - Vo(i-n));
            elseif (ismember(i,leboundin) || i==(n*m-2*n+2)) %leboundin, top-left internal node
                 dudx(i) = (1/dx)*(Uo(i) - Uo(i-1));
                 dudy(i) = (1/(2*dy))*(3*Uo(i) - 4*Uo(i-n) + Uo(i-2*n));
                 dvdx(i) = (1/dx)*(Vo(i) - Vo(i-1));
                 dvdy(i) = (1/(2*dy))*(3*Vo(i) - 4*Vo(i-n) + Vo(i-2*n));
            elseif i==(n+2) %bottom left internal node
                 dudx(i) = (1/dx)*(Uo(i) - Uo(i-1));
                 dudy(i) = (1/dy)*(Uo(i) - Uo(i-n));
                 dvdx(i) = (1/dx)*(Vo(i) - Vo(i-1));
                 dvdy(i) = (1/dy)*(Vo(i) - Vo(i-n));
            else %riboundin,upboundin,top-right internal node, pure internal nodes
                 dudx(i) = (1/(2*dx))*(3*Uo(i) - 4*Uo(i-1) + Uo(i-2));
                 dudy(i) = (1/(2*dy))*(3*Uo(i) - 4*Uo(i-n) + Uo(i-2*n));
                 dvdx(i) = (1/(2*dx))*(3*Vo(i) - 4*Vo(i-1) + Vo(i-2));
                 dvdy(i) = (1/(2*dy))*(3*Vo(i) - 4*Vo(i-n) + Vo(i-2*n));
            end
        elseif Uo(i) <= 0 && Vo(i) > 0  %which means forward upwind in x, and backward upwind in y
            %[doboundin, bottom-left internal node], second order in x, first order in y
            %[riboundin, top-right internal node], first order in x, second order in y
            %[bottom-right internal node], first order in both directions
            %[leboundin,upboundin,top-left internal node, pure internal nodes], second order in x and y
            if ismember(i,bound)
                continue
            elseif (ismember(i,doboundin) || i==(n+2)) %[doboundin, bottom-left internal node]
                dudx(i) = (1/(2*dx))*(-3*Uo(i) + 4*Uo(i+1) - Uo(i+2));
                dudy(i) = (1/dy)*(Uo(i) - Uo(i-n));
                dvdx(i) = (1/(2*dx))*(-3*Vo(i) + 4*Vo(i+1) - Vo(i+2));
                dvdy(i) = (1/dy)*(Vo(i) - Vo(i-n));
            elseif (ismember(i,riboundin) || i==(n*m-n-1)) %[riboundin, top-right internal node]
                dudx(i) = (1/dx)*(-Uo(i) + Uo(i+1));
                dudy(i) = (1/(2*dy))*(3*Uo(i) - 4*Uo(i-n) + Uo(i-2*n));
                dvdx(i) = (1/dx)*(-Vo(i) + Vo(i+1));
                dvdy(i) = (1/(2*dy))*(3*Vo(i) - 4*Vo(i-n) + Vo(i-2*n));
            elseif i==(2*n-1) %[bottom-right internal node]
                dudx(i) = (1/dx)*(-Uo(i) + Uo(i+1));
                dudy(i) = (1/dy)*(Uo(i) - Uo(i-n));
                dvdx(i) = (1/dx)*(-Vo(i) + Vo(i+1));
                dvdy(i) = (1/dy)*(Vo(i) - Vo(i-n));
            else %[leboundin,upboundin,top-left internal node, pure internal nodes]
                dudx(i) = (1/(2*dx))*(-3*Uo(i) + 4*Uo(i+1) - Uo(i+2));
                dudy(i) = (1/(2*dy))*(3*Uo(i) - 4*Uo(i-n) + Uo(i-2*n));
                dvdx(i) = (1/(2*dx))*(-3*Vo(i) + 4*Vo(i+1) - Vo(i+2));
                dvdy(i) = (1/(2*dy))*(3*Vo(i) - 4*Vo(i-n) + Vo(i-2*n));
            end
        elseif Uo(i) > 0 && Vo(i) <= 0  %which means backward upwind in x, and forward upwind in y
            %[upboundin, top-right internal node], second order in x, first order in y
            %[leboundin, bottom-left internal node], first order in x, second order in y
            %[top-left internal node], first order in both directions
            %[riboundin,doboundin,bottom-right internal node, pure internal nodes], second order in x and y
            if ismember(i,bound)
                continue
            elseif (ismember(i,upboundin) || i==(n*m-n-1)) %[upboundin, top-right internal node]
                dudx(i) = (1/(2*dx))*(3*Uo(i) - 4*Uo(i-1) + Uo(i-2));
                dudy(i) = (1/dy)*(-Uo(i) + Uo(i+n));
                dvdx(i) = (1/(2*dx))*(3*Vo(i) - 4*Vo(i-1) + Vo(i-2));
                dvdy(i) = (1/dy)*(-Vo(i) + Vo(i+n));
            elseif (ismember(i,leboundin) || i==(n+2)) %[leboundin, bottom-left internal node]
                dudx(i) = (1/dx)*(Uo(i) - Uo(i-1));
                dudy(i) = (1/(2*dy))*(-3*Uo(i) + 4*Uo(i+n) - Uo(i+2*n));
                dvdx(i) = (1/dx)*(Vo(i) - Vo(i-1));
                dvdy(i) = (1/(2*dy))*(-3*Vo(i) + 4*Vo(i+n) - Vo(i+2*n));
            elseif i==(n*m-2*n+2) %[top-left internal node]
                dudx(i) = (1/dx)*(Uo(i) - Uo(i-1));
                dudy(i) = (1/dy)*(-Uo(i) + Uo(i+n));
                dvdx(i) = (1/dx)*(Vo(i) - Vo(i-1));
                dvdy(i) = (1/dy)*(-Vo(i) + Vo(i+n));
            else %[riboundin,doboundin,bottom-right internal node, pure internal nodes]
                dudx(i) = (1/(2*dx))*(3*Uo(i) - 4*Uo(i-1) + Uo(i-2));
                dudy(i) = (1/(2*dy))*(-3*Uo(i) + 4*Uo(i+n) - Uo(i+2*n));
                dvdx(i) = (1/(2*dx))*(3*Vo(i) - 4*Vo(i-1) + Vo(i-2));
                dvdy(i) = (1/(2*dy))*(-3*Vo(i) + 4*Vo(i+n) - Vo(i+2*n));
            end
        else %which means forward upwind in both directions
            %[upboundin, top-left internal node], second order in x, first order in y
            %[riboundin, bottom-right internal node], first order in x, second order in y
            %[top-right internal node], first order in both directions
            %[leboundin,doboundin,bottom-left internal node, pure internal nodes], second order in x and y
            if ismember(i,bound)
                continue
            elseif (ismember(i,upboundin) || i==(n*m-2*n+2)) %[upboundin, top-left internal node]
                dudx(i) = (1/(2*dx))*(-3*Uo(i) + 4*Uo(i+1) - Uo(i+2));
                dudy(i) = (1/dy)*(-Uo(i) + Uo(i+n));
                dvdx(i) = (1/(2*dx))*(-3*Vo(i) + 4*Vo(i+1) - Vo(i+2));
                dvdy(i) = (1/dy)*(-Vo(i) + Vo(i+n));
            elseif (ismember(i,riboundin) || i==(2*n-1)) %[riboundin, bottom-right internal node]
                dudx(i) = (1/dx)*(-Uo(i) + Uo(i+1));
                dudy(i) = (1/(2*dy))*(-3*Uo(i) + 4*Uo(i+n) - Uo(i+2*n));
                dvdx(i) = (1/dx)*(-Vo(i) + Vo(i+1));
                dvdy(i) = (1/(2*dy))*(-3*Vo(i) + 4*Vo(i+n) - Vo(i+2*n));
            elseif i==(n*m-n-1) %[top-right internal node]
                dudx(i) = (1/dx)*(-Uo(i) + Uo(i+1));
                dudy(i) = (1/dy)*(-Uo(i) + Uo(i+n));
                dvdx(i) = (1/dx)*(-Vo(i) + Vo(i+1));
                dvdy(i) = (1/dy)*(-Vo(i) + Vo(i+n));
            else %[leboundin,doboundin,bottom-left internal node, pure internal nodes]
                dudx(i) = (1/(2*dx))*(-3*Uo(i) + 4*Uo(i+1) - Uo(i+2));
                dudy(i) = (1/(2*dy))*(-3*Uo(i) + 4*Uo(i+n) - Uo(i+2*n));
                dvdx(i) = (1/(2*dx))*(-3*Vo(i) + 4*Vo(i+1) - Vo(i+2));
                dvdy(i) = (1/(2*dy))*(-3*Vo(i) + 4*Vo(i+n) - Vo(i+2*n));
            end
        end
    end
            
end