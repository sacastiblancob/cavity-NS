% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Initial_condition_40_251_91_Re200.mat

%Space differentials
dx = (xmax-xmin)/(n-1);           %X diferential
dy = (ymax-ymin)/(m-1);           %Y diferential
x = xmin:dx:xmax;                %x vector
y = ymin:dy:ymax;                %y vector
X = kron(ones(1,m),x);
Y = kron(y,ones(1,n));
X = reshape(X,n,m)';
Y = reshape(Y,n,m)';
Y = Y(m:-1:1,:);

U = reshape(Uo,n,m);
U = U';
U = U(m:-1:1,:);
V = reshape(Vo,n,m);
V = V';
V = V(m:-1:1,:);

%Adding aditional coordinate to the right in U and V for posterior
%interpolation process
UI =[U(:,1) U zeros(m,1)];
VI =[zeros(m,1) V zeros(m,1)];
XI =[(xmin-2)*ones(m,1) X (xmax+2)*ones(m,1)];
YI =[Y(:,1) Y Y(:,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL DG MESH

%Number of number of subdomains in X and Y
Nxg = 20;
Nyg = 10;

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

% plot(XG,YG,'o','Color','k')
% axis equal

%CFL, NT
% CFL = 3;
% mindx = min(abs(XG(1,2:end)-XG(1,1:Nx*Nxg)));
% mindy = min(abs(YG(2:end,1)-YG(1:Ny*Nyg,1)));
% maxvel = max(max(max(abs(U),abs(V))));
% dt = CFL*(min(mindx,mindy)/maxvel);
% % dt = 0.05;
% NT = 35;

CFL = 3;
maxvel = max(max(max(abs(U),abs(V))));
dt = CFL*(min(dx,dy)/maxvel);
NT = 35;

% for j=1:Nyg
%     for i=1:Nxg
%         XDGgt0 = XG((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1);
%         YDGgt0 = YG((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1);
%         
%         % FTLE
%         % FTLE
%         Jx = JX(i);
%         Jy = JY(j);
%         [FTLEf,FTLEb,XDGgt2,YDGgt2,XDGgt1,YDGgt1,XDGgt0,YDGgt0] = ...
%             FTLE_dgelem(XI,YI,UI,VI,ep,et,XDGgt0,YDGgt0,dt,NT,Nx,Ny,Jx,Jy,XDD,YDD);
%         
%         %Storing results
%         XG1((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = XDGgt1;
%         YG1((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = YDGgt1;
%         XG2((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = XDGgt2;
%         YG2((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = YDGgt2;
%         FTLEF((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = FTLEf;
%         FTLEB((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = FTLEb;
% 
%     end
% end

for j=1:Nyg
    for i=1:Nxg
        XDGgt0 = XG((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1);
        YDGgt0 = YG((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1);
        
        % FTLE-forward
        Jx = JX(i);
        Jy = JY(j);
        [FTLEf,XDGgt1,YDGgt1] = ...
            FTLE_dgelem_f(XI,YI,UI,VI,XDGgt0,YDGgt0,dt,NT,Nx,Ny,Jx,Jy,XDD,YDD);

        %Storing results
        XG1f((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = XDGgt1;
        YG1f((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = YDGgt1;
        FTLEF((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = FTLEf;

        % FTLE-backward
        [FTLEb,XDGgt1,YDGgt1] = ...
            FTLE_dgelem_f(XI,YI,-UI,-VI,XDGgt0,YDGgt0,dt,NT,Nx,Ny,Jx,Jy,XDD,YDD);
        
        %Storing results
        XG1b((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = XDGgt1;
        YG1b((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = YDGgt1;        
        FTLEB((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = FTLEb;

    end
end

% load FTLE_3_35.mat
% save('FTLE_3_35.mat','XG','YG','XG1','YG1','FTLEF','FTLEB')

% figure(1)
% 
% % Plotting U velocity
% ax1 = subplot(2,1,1);
% surf(X,Y,U);
% view(0,90);
% % caxis([-0.4 0.8])
% colormap(ax1,jet)
% shading interp
% axis equal
% axis([xmin xmax ymin ymax -10 10]);
% colorbar
% title('U'); xlabel('X'); ylabel('Y');
% drawnow
% hold off
%     
% % Plotting V velocity
% ax2 = subplot(2,1,2);
% surf(X,Y,V);
% view(0,90);
% % caxis([0 4])
% colormap(ax2,jet)
% shading interp
% axis equal
% axis([xmin xmax ymin ymax -10 10]);
% colorbar
% title('V'); xlabel('X'); ylabel('Y');
% drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Computing LCS - Finite Differences

%CFL, NT
CFL = 3;
maxvel = max(max(max(abs(U),abs(V))));
dt = CFL*(min(dx,dy)/maxvel);
NT = 35;

XYD = zeros(m,n); YXD = zeros(m,n);
XD = X; YD=Y;
for t=1:NT
    [XD,YD] = AB3(XI,YI,UI,VI,XD,YD,dt);
end
XYD(2:m-1,2:n-1) = (1/(2*dy))*(XD(3:m,2:n-1) - XD(1:m-2,2:n-1));
YXD(2:m-1,2:n-1) = (1/(2*dx))*(YD(2:m-1,3:n) - YD(2:m-1,1:n-2));
XD(2:m-1,2:n-1) = (1/(2*dx))*(XD(2:m-1,3:n) - XD(2:m-1,1:n-2));
YD(2:m-1,2:n-1) = (1/(2*dy))*(YD(3:m,2:n-1) - YD(1:m-2,2:n-1));
TR = XD.^2 + YD.^2 + XYD.^2 + YXD.^2;
XD = (TR/2) + sqrt((TR/2).^2 - 1);
XD(:,1) = 0; XD(:,end) = 0; XD(1,:) = 0; XD(end,:) = 0;
XD = real(XD);
XDF = (1/dt)*log(sqrt(XD));
lim=35;

figure(2)

% Plotting FTLE-f-DG
ax3 = subplot(2,2,1);
surf(XG,YG,FTLEF);
% surf(XG1,YG1,FTLEF);
view(0,90);
caxis([-1 0.2])
% caxis([-2 30])
colormap(ax3,jet)
shading interp
axis equal
% axis([xmin xmax ymin ymax -10 10]);
colorbar
title('FTLE-f'); xlabel('X'); ylabel('Y');
drawnow

% Plotting FTLE-b-DG
ax4 = subplot(2,2,2);
surf(XG,YG,FTLEB);
% surf(XG1,YG1,FTLEB);
view(0,90);
caxis([-1 0.2])
% caxis([-2 30])
colormap(ax4,jet)
shading interp
axis equal
% axis([xmin xmax ymin ymax -10 10]);
colorbar
title('FTLE-b'); xlabel('X'); ylabel('Y');
drawnow

% Plotting FTLE-f
ax3 = subplot(2,2,3);
surf(X,Y,XDF);
view(0,90);
caxis([0 lim])
colormap(ax3,jet)
shading interp
axis equal
axis([xmin xmax ymin ymax -100 100]);
colorbar
title('FTLE-f'); xlabel('X'); ylabel('Y');
drawnow

%%%% BACKWARD
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
XD = real(XD);
XDF = (1/dt)*log(sqrt(XD));

% Plotting FTLE-b
ax4 = subplot(2,2,4);
surf(X,Y,XDF);
view(0,90);
caxis([0 lim])
colormap(ax4,jet)
shading interp
axis equal
axis([xmin xmax ymin ymax -100 100]);
colorbar
title('FTLE-b'); xlabel('X'); ylabel('Y');
drawnow














