function [ppgx,ppgy,ppgpx,ppgpy] = param(ep,et,xb1,xb2,xb3,xb4,yb1,yb2,yb3,yb4)
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
% ppgx, parametrizations for the 4 boundaries of x
% ppgx, parametrizations for the 4 boundaries of y
% ppgpx, derivatives of parametrizations for the 4 boundaries of x
% ppgpy, derivatives of parametrizations for the 4 boundaries of y
%
% Based on: DG-FTLE: Lagrangian coherent structures with high-order
%           discontinuous-Galerkin methods
%               Daniel A. Nelson, Gustaaf B. Jacobs, 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PARAMETRIZATION OF ELEMENT BOUNDARIES

%Gamma functions for boundaries in X and Y
ppg1x = spline(ep,xb1); ppg1y = spline(ep,yb1);
ppg2x = spline(et,xb2); ppg2y = spline(et,yb2);
ppg3x = spline(ep,xb3); ppg3y = spline(ep,yb3);
ppg4x = spline(et,xb4); ppg4y = spline(et,yb4);

%Derivatives of 
ppg1px = fnder(ppg1x,1); ppg1py = fnder(ppg1y,1);
ppg2px = fnder(ppg2x,1); ppg2py = fnder(ppg2y,1);
ppg3px = fnder(ppg3x,1); ppg3py = fnder(ppg3y,1);
ppg4px = fnder(ppg4x,1); ppg4py = fnder(ppg4y,1);

%Storing structures
ppgx=[ppg1x;ppg2x;ppg3x;ppg4x];
ppgy=[ppg1y;ppg2y;ppg3y;ppg4y];
ppgpx=[ppg1px;ppg2px;ppg3px;ppg4px];
ppgpy=[ppg1py;ppg2py;ppg3py;ppg4py];

end