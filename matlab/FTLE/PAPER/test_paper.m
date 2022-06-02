%%%% TEST
%
% Based on: DG-FTLE: Lagrangian coherent structures with high-order
%           discontinuous-Galerkin methods
%               Daniel A. Nelson, Gustaaf B. Jacobs, 2015


%Space differentials
xmax = 2; xmin = 0;
ymax = 1; ymin = 0;
% xmax = 1; xmin = -1;
% ymax = 1; ymin = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL DG MESH

%Number of number of subdomains in X and Y
Nxg = 32;
Nyg = 16;
% Nxg = 10;
% Nyg = 10;

%Lengths-of-elements vectors (here all have same length)
LXG = ((xmax-xmin)/Nxg)*ones(1,Nxg);
LYG = ((ymax-ymin)/Nyg)*ones(Nyg,1);

%Jacobians
JX = 1./LXG;
JY = 1./LYG;

%Local polynomial orders (all elements will have same discretization)
Nx = 12;
Ny = 12;

%Local-master element coordinates
ep = (0.5 + 0.5*JacobiGL(0,0,Nx))';
et = (0.5 + 0.5*JacobiGL(0,0,Ny));

% Derivative matrices
XDD = dmatrix(ep,Nx);
YDD = dmatrix(et,Ny);

%Coordinates of X-direction element boundaries
XBH = zeros(1,Nx*Nxg+1);

%Coordinates of Y-direction element boundaries
YBV = zeros(Ny*Nyg+1,1);

for i=1:Nxg
    XBH((i-1)*Nx+1:i*Nx+1) = LXG(i)*ep + (xmin + sum(LXG(1:i-1)));
end
for i=1:Nyg
    YBV((i-1)*Ny+1:i*Ny+1) = LYG(i)*et + (ymin + sum(LYG(1:i-1)));
end

%Global grid
[XG,YG] = meshgrid(XBH,YBV);
XG1 = XG; YG1 = YG;
XG2 = XG; YG2 = YG;
FTLEF = zeros(size(XG1)); FTLEB = zeros(size(XG1));

%CFL, NT
CFL = 3;
mindx = min(abs(XG(1,2:end)-XG(1,1:Nx*Nxg)));
mindy = min(abs(YG(2:end,1)-YG(1:Ny*Nyg,1)));
% t1 = pi/2 + pi/8;
t1 = 0;
[U,V] = f_vel(XG,YG,t1);
% streamslice(XG,YG,U,V)
A = 0.1;
maxvel = pi*A*sin(pi*0.5);
% maxvel = max(max(max(abs(U),abs(V))));
% dt = CFL*(min(mindx,mindy)/maxvel);
dt = .00125;
NT = 1600;

% % plot(XG,YG,'o','Color','k')
% % axis equal

for j=1:Nyg
    for i=1:Nxg
        XDGgt0 = XG((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1);
        YDGgt0 = YG((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1);
        
        % FTLE
        % FTLE
        Jx = JX(i);
        Jy = JY(j);
        [FTLEf,FTLEb,XDGgt2,YDGgt2,XDGgt1,YDGgt1] = ...
            FTLE_dgelem_fvel(ep,et,XDGgt0,YDGgt0,dt,NT,t1,Nx,Ny,XDD,YDD);
        
        %Storing results
        XG1((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = XDGgt1;
        YG1((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = YDGgt1;
        XG2((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = XDGgt2;
        YG2((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = YDGgt2;
        FTLEF((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = FTLEf;
        FTLEB((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = FTLEb;

    end
end
% load Results_gyre_3_20.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%5
% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)

% Plotting U velocity
ax1 = subplot(3,1,1);
streamslice(XG,YG,U,V);
view(0,90);
% caxis([-0.4 0.8])
axis equal
axis([xmin xmax ymin ymax -100 100]);
title('Stream'); xlabel('X'); ylabel('Y');
drawnow
hold off
    
% Plotting FTLE-f
ax3 = subplot(3,1,2);
surf(XG1,YG1,FTLEF);
view(0,90);
% caxis([-0.4 0.4])
shading interp
colormap jet
axis equal
axis([xmin xmax ymin ymax min(min(FTLEF)) max(max(FTLEF))]);
colorbar
title('FTLE-f'); xlabel('X'); ylabel('Y');
drawnow

% Plotting FTLE-b
ax4 = subplot(3,1,3);
surf(XG1,YG1,FTLEB);
view(0,90);
% caxis([-0.4 0.4])
shading interp
colormap jet
axis equal
axis([xmin xmax ymin ymax min(min(FTLEB)) max(max(FTLEB))]);
colorbar
title('FTLE-b'); xlabel('X'); ylabel('Y');
drawnow








