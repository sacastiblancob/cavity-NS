function [XDGg3,YDGg3] = AB3(X,Y,U,V,XDGg,YDGg,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the particle tracing with Adam-Bashfort-3 for a
% time step dt of location initialized at XDGg and YDGg over velocity field
% U,V.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = dt/3;

%Euler
% Velocity field interpolated at particle locations
UDG0 = interp2(X,Y,U,XDGg,YDGg,'makima');
VDG0 = interp2(X,Y,V,XDGg,YDGg,'makima');

% The particles are advected one step in time
XDGg1 = XDGg + h*UDG0;
YDGg1 = YDGg + h*VDG0;

%AB-2
% Velocity field interpolated at particle locations
UDG1 = interp2(X,Y,U,XDGg1,YDGg1,'makima');
VDG1 = interp2(X,Y,V,XDGg1,YDGg1,'makima');

% The particles are advected one step in time
XDGg2 = XDGg1 + h*((3/2)*UDG1 - 0.5*UDG0);
YDGg2 = YDGg1 + h*((3/2)*VDG1 - 0.5*VDG0);

%AB-3
% Velocity field interpolated at particle locations
UDG2 = interp2(X,Y,U,XDGg2,YDGg2,'makima');
VDG2 = interp2(X,Y,V,XDGg2,YDGg2,'makima');

% The particles are advected one step in time
XDGg3 = XDGg2 + h*((23/12)*UDG2 - (16/12)*UDG1 + (5/12)*UDG0);
YDGg3 = YDGg2 + h*((23/12)*VDG2 - (16/12)*VDG1 + (5/12)*VDG0);


end