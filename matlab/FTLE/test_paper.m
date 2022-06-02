%%%% TEST
%
% Based on: DG-FTLE: Lagrangian coherent structures with high-order
%           discontinuous-Galerkin methods
%               Daniel A. Nelson, Gustaaf B. Jacobs, 2015


%Space differentials
xmax = 1; xmin = 0;
ymax = 1; ymin = 0;

n = 51;
m = 51;

dx = (xmax-xmin)/(n-1);           %X diferential
dy = (ymax-ymin)/(m-1);           %Y diferential
x = xmin:dx:xmax;                %x vector
y = ymin:dy:ymax;                %y vector
X = kron(ones(1,m),x);
Y = kron(y,ones(1,n));
X = reshape(X,n,m)';
Y = reshape(Y,n,m)';
% Y = Y(m:-1:1,:);

A = 1;
U = -pi*A.*sin(pi*X).*cos(pi*Y);
V = pi*A.*sin(pi*Y).*cos(pi*X);
% 
% streamslice(X,Y,U,V)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Computing LCS - DG

Nx = 24;
Ny = 24;
ep = (0.5 + 0.5*JacobiGL(0,0,Nx))';
et = (0.5 + 0.5*JacobiGL(0,0,Ny));

x0dg = 0;
y0dg = 0;
lxdg = 1;
lydg = 1;
Jx = 1/lxdg;
Jy = 1/lydg;

XDGl = kron(ones(Ny+1,1),ep);
YDGl = kron(et,ones(1,Nx+1));

XDGg = lxdg*XDGl + x0dg;
YDGg = lydg*YDGl + y0dg;
XDGgt0 = XDGg;
YDGgt0 = YDGg;

% Derivative matrices
XDD = dmatrix(ep,Nx);
YDD = dmatrix(et,Ny);

%CFL, NT
CFL = 1;
mindx = min(abs(X(1,2:end)-X(1,1:n-1)));
mindy = min(abs(Y(2:end,1)-Y(1:m-1,1)));
maxvel = max(max(max(abs(U),abs(V))));
dt = CFL*(min(mindx,mindy)/maxvel);
NT = 10;

% FTLE
[FTLEf,FTLEb,XDGgt2,YDGgt2,XDGgt1,YDGgt1,XDGgt0,YDGgt0] = ...
    FTLE_dgelem(X,Y,U,V,ep,et,XDGgt0,YDGgt0,dt,NT,Nx,Ny,Jx,Jy,XDD,YDD);
% [FTLEf,FTLEb,XDGgt2,YDGgt2,XDGgt1,YDGgt1,XDGgt0,YDGgt0] = ...
%     FTLE_dgelem2(X,Y,U,V,ep,et,XDGgt0,YDGgt0,dt,NT,Nx,Ny,Jx,Jy,XDD,YDD);
% [FTLEf,FTLEb,XDGgt2,YDGgt2,XDGgt1,YDGgt1,XDGgt0,YDGgt0] = ...
%     FTLE_dgelem3(X,Y,U,V,ep,et,XDGgt0,YDGgt0,dt,NT,Nx,Ny,Jx,Jy,XDD,YDD);

% Plotting FTLE-f
figure(1)
surf(XDGgt1,YDGgt1,FTLEf);
view(0,90);
shading interp
axis equal
% axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -1000 1000]);
colorbar
title('FTLE-f'); xlabel('X'); ylabel('Y');
drawnow

% Plotting FTLE-b
figure(2)
surf(XDGgt1,YDGgt1,FTLEb);
view(0,90);
shading interp
axis equal
% axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -1000 1000]);
colorbar
title('FTLE-b'); xlabel('X'); ylabel('Y');
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Computing LCS - Finite Differences

%CFL, NT
CFL = 1;
mindx = min(abs(X(1,2:end)-X(1,1:n-1)));
mindy = min(abs(Y(2:end,1)-Y(1:m-1,1)));
maxvel = max(max(max(abs(U),abs(V))));
dt = CFL*(min(mindx,mindy)/maxvel);
NT = 10;

XYD = zeros(m,n); YXD = zeros(m,n);
XD = X; YD=Y;
for t=1:NT
    [XD,YD] = AB3(X,Y,U,V,XD,YD,dt);
end
XYD(2:m-1,2:n-1) = (1/(2*dy))*(XD(3:m,2:n-1) - XD(1:m-2,2:n-1));
YXD(2:m-1,2:n-1) = (1/(2*dx))*(YD(2:m-1,3:n) - YD(2:m-1,1:n-2));
XD(2:m-1,2:n-1) = (1/(2*dx))*(XD(2:m-1,3:n) - XD(2:m-1,1:n-2));
YD(2:m-1,2:n-1) = (1/(2*dy))*(YD(3:m,2:n-1) - YD(1:m-2,2:n-1));
TR = XD.^2 + YD.^2 + XYD.^2 + YXD.^2;
XD = (TR/2) + sqrt((TR/2).^2 - 1);
XD(:,1) = 0; XD(:,end) = 0; XD(1,:) = 0; XD(end,:) = 0; 
XDF = (1/dt)*log(sqrt(XD));

% Plotting FTLE-f
figure(3)
surf(X,Y,XDF);
view(0,90);
shading interp
axis equal
% axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -1000 1000]);
colorbar
title('FTLE-f'); xlabel('X'); ylabel('Y');
drawnow

%CFL, NT
CFL = 1;
mindx = min(abs(X(1,2:end)-X(1,1:n-1)));
mindy = min(abs(Y(2:end,1)-Y(1:m-1,1)));
maxvel = max(max(max(abs(U),abs(V))));
dt = CFL*(min(mindx,mindy)/maxvel);
NT = 10;

XYD = zeros(m,n); YXD = zeros(m,n);
XD = X; YD=Y;
for t=1:NT
    [XD,YD] = AB3(X,Y,-U,-V,XD,YD,dt);
end
XYD(2:m-1,2:n-1) = (1/(2*dy))*(XD(3:m,2:n-1) - XD(1:m-2,2:n-1));
YXD(2:m-1,2:n-1) = (1/(2*dx))*(YD(2:m-1,3:n) - YD(2:m-1,1:n-2));
XD(2:m-1,2:n-1) = (1/(2*dx))*(XD(2:m-1,3:n) - XD(2:m-1,1:n-2));
YD(2:m-1,2:n-1) = (1/(2*dy))*(YD(3:m,2:n-1) - YD(1:m-2,2:n-1));
TR = XD.^2 + YD.^2 + XYD.^2 + YXD.^2;
XD = (TR/2) + sqrt((TR/2).^2 - 1);
XD(:,1) = 0; XD(:,end) = 0; XD(1,:) = 0; XD(end,:) = 0; 
XDB = (1/dt)*log(sqrt(XD));

% Plotting FTLE-b
figure(4)
surf(X,Y,XDB);
view(0,90);
shading interp
axis equal
% axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -1000 1000]);
colorbar
title('FTLE-b'); xlabel('X'); ylabel('Y');
drawnow




