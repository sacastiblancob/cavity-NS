function [dudx,dvdy] = grad2D(Uo,Vo,dx,dy,n,m,upbound,dobound,...
    lebound,ribound)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the terms du/dx, dv/dy with O(h^2) centered
% differences (if posible)
%
% Poisson equation
%
%   dP2/d2x + dP2/d2y = -(rho/dt)*(du*/dx + dv*/dy)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entries
%   Uo, vector with velocities in X direction
%   Uv, vector with velocities in Y direction
%   dx, differential in X
%   dy, differential in Y
%   n, number of nodes in X direction
%   m, number of nodes in Y direction
%   upbound, indices of the top boundary (ub)
%   dobound, indices of the bottom boundary (db)
%   lebound, indices of the left boundary (lb)
%   ribound, indices of the right boundary (rb)
%      *[upbound, dobound, lebound, ribound] € bound
%   boundin, internal boundaries
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
% upbound(ub), dobound(db), lebound(lb), ribound(rb), boundint(bi)
% internal nodes(i)
%
%    ulb  ub   ub   ub   ub   urb
%
%    lb   i    i    i    i    rb
%
%    lb   i    i    i    i    rb
%
%    lb   i    i    i    i    rb
%
%    dlb  db   db   db   db   drb
%
% dul --> top-left node
% ddl --> bottom-left node
% ddr --> bottom-right node
% dur --> top-right node
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %Allocating memory for the result
    dudx = zeros(n*m,1);
    dvdy = zeros(n*m,1);
    
    %Loop
    for i=1:n*m
        %
        % Derivatives in the corners
        %
        if i==1
            dudx(i) = (1/dx)*(Uo(i+1) - Uo(i));
            dvdy(i) = (1/dy)*(Vo(i+n) - Vo(i));
        elseif i == n
            dudx(i) = (1/dx)*(Uo(i) - Uo(i-1));
            dvdy(i) = (1/dy)*(Vo(2*i) - Vo(i));
        elseif i == n*m-n+1
            dudx(n*m-n+1) = (1/dx)*(Uo(n*m-n+2) - Uo(n*m-n+1));
            dvdy(n*m-n+1) = (1/dy)*(Vo(n*m-n+1) - Vo(n*m-2*n+1));
        elseif i == n*m
            dudx(n*m) = (1/dx)*(Uo(n*m) - Uo(n*m-1));
            dvdy(n*m) = (1/dy)*(Vo(n*m) - Vo(n*m-n));
        %
        % Derivatives in the boundaries
        %
        elseif ismember(i,upbound)
            dudx(i) = (1/(2*dx))*(Uo(i+1) - Uo(i-1));
            dvdy(i) = (1/dy)*(Vo(i) - Vo(i-n));
        elseif ismember(i,dobound)
            dudx(i) = (1/(2*dx))*(Uo(i+1) - Uo(i-1));
            dvdy(i) = (1/dy)*(Vo(i+n) - Vo(i));
        elseif ismember(i,lebound)
            dudx(i) = (1/dx)*(Uo(i+1) - Uo(i));
            dvdy(i) = (1/(2*dy))*(Vo(i+n) - Vo(i-n));
        elseif ismember(i,ribound)
            dudx(i) = (1/dx)*(Uo(i) - Uo(i-1));
            dvdy(i) = (1/(2*dy))*(Vo(i+n) - Vo(i-n));
        %
        % Derivatives in the internal portion
        %
        else
            dudx(i) = (1/(2*dx))*(Uo(i+1) - Uo(i-1));
            dvdy(i) = (1/(2*dy))*(Vo(i+n) - Vo(i-n));
        end
    end

end