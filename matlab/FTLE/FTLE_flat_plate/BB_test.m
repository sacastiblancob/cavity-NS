% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Space differentials
xmin = -0.5;
ymin = 0;
xmax = 2;
ymax = 0.25;

n = 251;
m = 151;
x = xmin:(xmax-xmin)/(n-1):xmax;
y = ymin:(ymax-ymin)/(m-1):ymax;
[X,Y] = meshgrid(x,y);

[U,V] = f_vel(X,Y,0);
% streamslice(X,Y,U,V)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % TRANSPORT OF MATERIAL LINES CLOSE TO THE WALL
% yw = y(1:5:5*6);
% [XM,YM] = meshgrid(x,yw);
% 
% % % Plotting Material Lines
% % ax2 = subplot(3,1,3);
% % streamslice(XS,YS,US,VS);
% % hold on
% % plot(XM,YM,'.','Color','k');
% % axis equal
% % axis([xmins xmaxs ymins ymaxs -10 10]);
% % title('Material lines'); xlabel('X'); ylabel('Y');
% % drawnow
% 
% NT=200;
% dtp = 0.01;
% for t=1:NT
%     [XM,YM] = AB3_fvel(XM,YM,dtp,0);
% 
% %     % Plotting Material Lines
% %     ax2 = subplot(3,1,3);
% %     plot(XM,YM,'.','Color','k');
% %     axis equal
% %     axis([xmins xmaxs ymins ymaxs -10 10]);
% %     title('Material lines'); xlabel('X'); ylabel('Y');
% %     drawnow
% %     pause
% end
% 
% % Plotting Material Lines
% streamslice(X,Y,U,V);
% hold on
% plot(XM,YM,'.','Color','k');
% axis equal
% axis([xmin xmax ymin ymax -10 10]);
% title('Material lines'); xlabel('X'); ylabel('Y');
% drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL DG MESH

%Number of number of subdomains in X and Y
Nxg = 20;
Nyg = 1;

ymaxs = 0.04;
xmins = 0;
%Lengths-of-elements vectors (here all have same length)
LXG = ((xmax-xmins)/Nxg)*ones(1,Nxg);
LYG = ((ymaxs-ymin)/Nyg)*ones(Nyg,1);

%Jacobians
JX = 1./LXG;
JY = 1./LYG;

%Local polynomial orders (all elements will have same discretization)
Nx = 14;
Ny = 14;

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
    XBH((i-1)*Nx+1:i*Nx+1) = LXG(i)*ep + (xmins + sum(LXG(1:i-1)));
end
for i=1:Nyg
    YBV((i-1)*Ny+1:i*Ny+1) = LYG(i)*et + (ymin + sum(LYG(1:i-1)));
end

%Global grid
[XG,YG] = meshgrid(XBH,YBV);
XG1 = XG; YG1 = YG;
XG2 = XG; YG2 = YG;
% FTLEF = zeros(size(XG1));
FTLEB = zeros(size(XG1));
K0 = zeros(size(XG1));
DXDEP = zeros(size(XG1));
DYDEP = zeros(size(XG1));
DXDEP2 = zeros(size(XG1));
DYDEP2 = zeros(size(XG1));

% plot(XG,YG,'o','Color','k')
% axis equal

%CFL, NT
% CFL = 3;
% mindx = min(abs(XG(1,2:end)-XG(1,1:Nx*Nxg)));
% mindy = min(abs(YG(2:end,1)-YG(1:Ny*Nyg,1)));
% maxvel = max(max(max(abs(U),abs(V))));
% dt = CFL*(min(mindx,mindy)/maxvel);
dt = 0.005;
NT = 0;
NT1 = 400;

for j=1:Nyg
    for i=1:Nxg
        XDGgt0 = XG((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1);
        YDGgt0 = YG((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1);
        
        % FTLE-forward
        Jx = JX(i);
        Jy = JY(j);
%         [FTLEf,XDGgt1,YDGgt1] = ...
%             FTLE_dgelem_f(X,Y,U,V,XDGgt0,YDGgt0,dt,NT,Nx,Ny,Jx,Jy,XDD,YDD);
% 
%         %Storing results
%         XG1f((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = XDGgt1;
%         YG1f((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = YDGgt1;
%         FTLEF((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = FTLEf;

%         % FTLE-backward
%         [FTLEb,XDGgt1,YDGgt1] = ...
%             FTLE_dgelem_f(X,Y,-U,-V,XDGgt0,YDGgt0,dt,NT,Nx,Ny,Jx,Jy,XDD,YDD);
%         
%         %Storing results
%         XG1b((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = XDGgt1;
%         YG1b((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = YDGgt1;        
%         FTLEB((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = FTLEb;

%         % Curvature scalars
%         [k0l,dxdep,dydep,XDGgt1,YDGgt1,dxdep2,dydep2] = ...
%             curvk_dgelem(X,Y,U,V,ep,et,XDGgt0,YDGgt0,dt,NT,Nx,Ny,XDD,YDD);
%         
%         %Storing results
%         XG1b((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = XDGgt1;
%         YG1b((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = YDGgt1;
%         K0((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = k0l;
%         DXDEP((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = dxdep;
%         DYDEP((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = dydep;
%         DXDEP2((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = dxdep2;
%         DYDEP2((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = dydep2;

%         % Curvature scalars evolution for lagrangian Backbone
%         [Kt0,k0l,dxdep,dydep,XDGgt1,YDGgt1] = ...
%             lagcurvchange_dgelem_fvel(ep,et,XDGgt0,YDGgt0,dt,NT,NT1,Nx,Ny,XDD,YDD);
%         
%         %Storing results
%         XG1b((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = XDGgt1;
%         YG1b((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = YDGgt1;
%         K0((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = k0l;
%         KT0((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = Kt0;
%         DXDEP((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = dxdep;
%         DYDEP((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = dydep;

        % Curvature scalars evolution for eulerian Backbone
        [Kt0,k0l,dxdep,dydep,XDGgt1,YDGgt1] = ...
            eulcurvchange_dgelem_fvel(ep,et,XDGgt0,YDGgt0,dt,NT,NT1,Nx,Ny,XDD,YDD);
        
        %Storing results
        XG1b((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = XDGgt1;
        YG1b((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = YDGgt1;
        K0((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = k0l;
        KT0((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = Kt0;
        DXDEP((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = dxdep;
        DYDEP((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = dydep;

    end
end

%Extracting ridge from Kt0
ch = zeros(Nx*Nxg+1,1);
ch(1:Nx:Nx*Nxg+1) = 1;
ch = -ch+1;
ch = logical(ch);
xrid = zeros(Ny*Nyg+1,1);
yrid = zeros(Ny*Nyg+1,1);
zrid = zeros(Ny*Nyg+1,1);
for j=1:Ny*Nyg+1
    zrid(j) = max(KT0(j,ch));
    xrid(j) = XG1b(j,KT0(j,ch)==zrid(j));
    yrid(j) = YG1b(j,KT0(j,ch)==zrid(j));
end

figure(1)

% % Plotting U velocity
% ax1 = subplot(2,2,1);
% surf(X,Y,U);
% view(0,90);
% % caxis([-0.4 0.8])
% colormap(ax1,jet)
% shading interp
% % axis equal
% axis([xmin xmax ymin ymax -10 10]);
% colorbar
% title('U'); xlabel('X'); ylabel('Y');
% drawnow
% hold off
%     
% % Plotting V velocity
% ax2 = subplot(2,2,2);
% surf(X,Y,V);
% view(0,90);
% % caxis([0 4])
% colormap(ax2,jet)
% shading interp
% % axis equal
% axis([xmin xmax ymin ymax -10 10]);
% colorbar
% title('V'); xlabel('X'); ylabel('Y');
% drawnow

% % Plotting FTLE-b-DG
% ax3 = subplot(3,1,3);
% surf(XG,YG,FTLEB);
% hold on
% plot3(XM,YM,100*ones(size(XM))'.','Color','k');
% view(0,90);
% caxis([-2 -1])
% colormap(ax3,jet)
% shading interp
% axis equal
% axis([xmins xmaxs ymins ymaxs -100 100]);
% colorbar
% title('FTLE-b'); xlabel('X'); ylabel('Y');
% drawnow

% % Plotting Curvature Scalars
% ax3 = subplot(2,2,3);
ax3 = subplot(2,1,1);
quiver(XG1b,YG1b,DXDEP,DYDEP,'k')
hold on
plot(xrid,yrid,'r')
streamslice(X,Y,U,V,'noarrows')
axis([xmin xmax ymin ymax -1000 1000]);
% quiver(XG1b,YG1b,DXDEP2,DYDEP2,'k')
% axis equal

% surf(XG1b,YG1b,K0);
% hold on
% % quiver(XG1b,YG1b,DXDEP,DYDEP)
% % quiver(XG1b,YG1b,DXDEP2,DYDEP2,'k')
% % plot3(XM,YM,100*ones(size(XM))'.','Color','k');
% % plot3(XG1b,YG1b,0.05*ones(size(XG1b)),'.','Color','k');
% view(0,90);
% % caxis([-2 -1])
% colormap(ax3,jet)
% % colormap jet
% shading interp
% axis equal
% % axis([xmins xmaxs ymins ymaxs -100 100]);
% colorbar
% title('Curvature Scalars'); xlabel('X'); ylabel('Y');
% drawnow

% Plotting Curvature Scalars Evolution
% ax4 = subplot(2,2,4);
ax4 = subplot(2,1,2);
surf(XG1b,YG1b,KT0);
hold on
plot3(xrid,yrid,zrid,'r')
% quiver(XG1b,YG1b,DXDEP,DYDEP)
% quiver(XG1b,YG1b,DXDEP2,DYDEP2,'k')
% plot3(XM,YM,100*ones(size(XM))'.','Color','k');
% plot3(XG1b,YG1b,0.05*ones(size(XG1b)),'.','Color','k');
view(0,90);
% caxis([-2 -1])
colormap(ax4,jet)
% colormap jet
shading interp
% axis equal
axis([xmin xmax ymin ymax min(min(KT0)) max(max(KT0))]);
colorbar
title('Curvature Scalars Evo.'); xlabel('X'); ylabel('Y');
drawnow

% figure(2)
% quiver(XG1b,YG1b,DXDEP,DYDEP)
% hold on
% % quiver(XG1b,YG1b,DXDEP2,DYDEP2,'k')
% % axis equal







