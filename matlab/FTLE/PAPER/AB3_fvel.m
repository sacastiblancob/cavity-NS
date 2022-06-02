function [XDGg3,YDGg3] = AB3_fvel(XDGg,YDGg,dt,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the particle tracing with Adam-Bashfort-3 for a
% time step dt of location initialized at XDGg and YDGg over velocity field
% U,V.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = dt/3;

%Euler
% Velocity field interpolated at particle locations
[UDG0,VDG0] = f_vel(XDGg,YDGg,t);

% The particles are advected one step in time
XDGg1 = XDGg + h*UDG0;
YDGg1 = YDGg + h*VDG0;

%AB-2
% Velocity field interpolated at particle locations
[UDG1,VDG1] = f_vel(XDGg1,YDGg1,t+h);

% The particles are advected one step in time
XDGg2 = XDGg1 + h*((3/2)*UDG1 - 0.5*UDG0);
YDGg2 = YDGg1 + h*((3/2)*VDG1 - 0.5*VDG0);

%AB-3
% Velocity field interpolated at particle locations
[UDG2,VDG2] = f_vel(XDGg2,YDGg2,t+2*h);

% The particles are advected one step in time
XDGg3 = XDGg2 + h*((23/12)*UDG2 - (16/12)*UDG1 + (5/12)*UDG0);
YDGg3 = YDGg2 + h*((23/12)*VDG2 - (16/12)*VDG1 + (5/12)*VDG0);


end