
% load Initial_condition_60_251_91_Re200.mat
% 
% %Space differentials
% dx = (xmax-xmin)/(n-1);           %X diferential
% dy = (ymax-ymin)/(m-1);           %Y diferential
% x = xmin:dx:xmax;                %x vector
% y = ymin:dy:ymax;                %y vector
% X = kron(ones(1,m),x);
% Y = kron(y,ones(1,n));
% X = reshape(X,n,m)';
% Y = reshape(Y,n,m)';
% Y = Y(m:-1:1,:);
% 
% Uplot = reshape(Uo,n,m);
% Uplot = Uplot';
% Uplot = Uplot(m:-1:1,:);
% Vplot = reshape(Vo,n,m);
% Vplot = Vplot';
% Vplot = Vplot(m:-1:1,:);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % DG Element
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nx = 4;
% Ny = 4;
% ep = (0.5 + 0.5*JacobiGL(0,0,Nx))';
% et = (0.5 + 0.5*JacobiGL(0,0,Ny));
% 
% x0dg = 1.0;
% y0dg = 0.65;
% lxdg = 0.3;
% lydg = 0.3;
% 
% XDGl = kron(ones(Ny+1,1),ep);
% YDGl = kron(et,ones(1,Nx+1));
% 
% XDGg = lxdg*XDGl + x0dg;
% YDGg = lydg*YDGl + y0dg;
% 
% D = dmatrix(ep,Nx);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % STEP A - Particle tracing algorithm
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %3 Velocity field interpolated at particle locations
% UDG = interp2(X,Y,Uplot,XDGg,YDGg,'makima');
% VDG = interp2(X,Y,Vplot,XDGg,YDGg,'makima');
% 
% %4. The particles are advected one step in time
% XDGg1 = XDGg + dt*UDG;
% YDGg1 = YDGg + dt*VDG;
% 
% %parameter of parametrization is t
% 
% ppg1x = spline(ep,XDGg1(1,:)); ppg1y = spline(ep,YDGg1(1,:));
% ppg2x = spline(et,XDGg1(:,end)); ppg2y = spline(et,YDGg1(:,end));
% ppg3x = spline(ep,XDGg1(end,:)); ppg3y = spline(ep,YDGg1(end,:));
% ppg4x = spline(et,XDGg1(:,1)); ppg4y = spline(et,YDGg1(:,1));
% 
% ppg1px = fnder(ppg1x,1); ppg1py = fnder(ppg1y,1);
% ppg2px = fnder(ppg2x,1); ppg2py = fnder(ppg2y,1);
% ppg3px = fnder(ppg3x,1); ppg3py = fnder(ppg3y,1);
% ppg4px = fnder(ppg4x,1); ppg4py = fnder(ppg4y,1);
% 
% ep1 = 1;
% et1 = 1;
% g1 = [ppval(ppg1x,ep1); ppval(ppg1y,ep1)];
% g2 = [ppval(ppg2x,et1); ppval(ppg2y,et1)];
% g3 = [ppval(ppg3x,ep1); ppval(ppg3y,ep1)];
% g4 = [ppval(ppg4x,et1); ppval(ppg4y,et1)];
% x1 = [XDGg1(1,1); YDGg1(1,1)];
% x2 = [XDGg1(1,end); YDGg1(1,end)];
% x3 = [XDGg1(end,end); YDGg1(end,end)];
% x4 = [XDGg1(end,1); YDGg1(end,1)];
% xy1 = (1-et1)*g1 + et1*g3 + (1-ep1)*g4 + ep1*g2 - x1*(1-ep1)*(1-et1) - ...
%     x2*ep1*(1-et1) - x3*ep1*et1 - x4*(1-ep1)*et1;
% 
% %%%%%%%%%%%%%%%%%%
% % JACOBIAN
% %%%%
% g1p = [ppval(ppg1px,ep1); ppval(ppg1py,ep1)];
% g2p = [ppval(ppg2px,et1); ppval(ppg2py,et1)];
% g3p = [ppval(ppg3px,ep1); ppval(ppg3py,ep1)];
% g4p = [ppval(ppg4px,et1); ppval(ppg4py,et1)];
% dxdep = x1*(1-et1) - x2*(1-et1) - x3*et1 + x4*et1 + g1p*(1-et1) + g3p*et1 + g2 - g4;
% dxdet = x1*(1-ep1) + x2*ep1 - x3*ep1 - x4*(1-ep1) + g2p*ep1 + g4p*(1-ep1) - g1 + g3;
% 
% J = [dxdep dxdet];
% 
% xy0 = [XDGg1(5,6); YDGg1(5,6)];
% 
% F1 = (1-et1)*g1 + et1*g3 + (1-ep1)*g4 + ep1*g2 - x1*(1-ep1)*(1-et1) - ...
%     x2*ep1*(1-et1) - x3*ep1*et1 - x4*(1-ep1)*et1;
% 
% epet0 = [ep1;et1];
% 
% epet1 = epet0 - J\(F1-xy0);
% 
% % XDGg3 = zeros(Nx+1,Nx+1);
% % YDGg3 = zeros(Nx+1,Nx+1);
% % for j = 1:Nx+1
% %     for i = 1:Nx+1
% %         ep0 = ep(i);
% %         et0 = et(j);
% %         g1 = [ppval(ppgx(1),ep0); ppval(ppgy(1),ep0)];
% %         g2 = [ppval(ppgx(2),et0); ppval(ppgy(2),et0)];
% %         g3 = [ppval(ppgx(3),ep0); ppval(ppgy(3),ep0)];
% %         g4 = [ppval(ppgx(4),et0); ppval(ppgy(4),et0)];
% %         F1 = (1-et0)*g1 + et0*g3 + (1-ep0)*g4 + ep0*g2 - x1*(1-ep0)*(1-et0) - ...
% %             x2*ep0*(1-et0) - x3*ep0*et0 - x4*(1-ep0)*et0;
% %         XDGg3(i,j) = F1(1);
% %         YDGg3(i,j) = F1(2);
% %     end
% % end
% 
% 
% 
% % 
% % % xp = XDGg1(1,1):0.001:XDGg1(1,end);
% % % yp = ppval(ppg1,xp);
% % % ypp = ppval(ppg1p,xp) + mean(yp);
% % 
% figure(2)
% % Plotting U velocity
% % surf(X,Y,Uplot);
% hold on
% plot3(XDGgt0,YDGgt0,10*ones(size(XDGg)),'*','Color','b')
% plot3(XDGgt1,YDGgt1,10*ones(size(XDGg)),'*','Color','k')
% plot3(XDGgt2,YDGgt2,10*ones(size(XDGg)),'*','Color','g')
% % plot3(XDGgt3,YDGgt3,10*ones(size(XDGg)),'*','Color','r')
% % % plot3(XDGl,YDGl,10*ones(size(XDGg)),'*','Color','b')
% % % hold on
% % % plot3(XDGl1,YDGl1,10*ones(size(XDGg)),'*','Color','k')
% % plot3(xy1(1),xy1(2),10,'o','Color','k')
% view(0,90);
% caxis([-0.4 0.8])
% colormap jet
% shading interp
% axis equal
% axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -10 10]);
% colorbar
% title('U and DG_el');
% drawnow
% hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Initial_condition_60_251_91_Re200.mat

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % DG Element
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nx = 10;
% Ny = 10;
% ep = (0.5 + 0.5*JacobiGL(0,0,Nx))';
% et = (0.5 + 0.5*JacobiGL(0,0,Ny));
% 
% x0dg = 1.0;
% y0dg = 0.65;
% lxdg = 0.3;
% lydg = 0.3;
% Jx = 1/lxdg;
% Jy = 1/lydg;
% 
% XDGl = kron(ones(Ny+1,1),ep);
% YDGl = kron(et,ones(1,Nx+1));
% 
% XDGg = lxdg*XDGl + x0dg;
% YDGg = lydg*YDGl + y0dg;
% 
% D = dmatrix(ep,Nx);
% 
% CFL = 0.5; NT=6;
% 
% [FTLEf,FTLEb,XDGgt2,YDGgt2,XDGgt0,YDGgt0] = FTLE_dgelem(X,Y,U,V,ep,et,XDGg,YDGg,CFL,NT,Nx,Ny,Jx,Jy,D,D);
% 
% figure(2)
% 
% % Plotting U velocity
% ax1 = subplot(2,2,1);
% surf(X,Y,U);
% view(0,90);
% % caxis([-0.4 0.8])
% colormap(ax1,jet)
% shading interp
% axis equal
% axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -10 10]);
% colorbar
% title('U'); xlabel('X'); ylabel('Y');
% drawnow
% hold off
%     
% % Plotting V velocity
% ax2 = subplot(2,2,2);
% surf(X,Y,V);
% view(0,90);
% % caxis([-0.4 0.4])
% colormap(ax2,jet)
% shading interp
% axis equal
% axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -10 10]);
% colorbar
% title('V'); xlabel('X'); ylabel('Y');
% drawnow
% 
% % Plotting FTLE-f
% ax3 = subplot(2,2,3);
% surf(XDGgt2,YDGgt2,FTLEf);
% view(0,90);
% % caxis([-0.4 0.4])
% shading interp
% axis equal
% axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -1000 1000]);
% colorbar
% title('FTLE-f'); xlabel('X'); ylabel('Y');
% drawnow
% 
% % Plotting FTLE-b
% ax4 = subplot(2,2,4);
% surf(XDGgt0,YDGgt0,FTLEb);
% view(0,90);
% % caxis([-0.4 0.4])
% shading interp
% axis equal
% axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -1000 1000]);
% colorbar
% title('FTLE-b'); xlabel('X'); ylabel('Y');
% drawnow

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% load Initial_condition_60_251_91_Re200.mat
% 
% %Space differentials
% dx = (xmax-xmin)/(n-1);           %X diferential
% dy = (ymax-ymin)/(m-1);           %Y diferential
% x = xmin:dx:xmax;                %x vector
% y = ymin:dy:ymax;                %y vector
% X = kron(ones(1,m),x);
% Y = kron(y,ones(1,n));
% X = reshape(X,n,m)';
% Y = reshape(Y,n,m)';
% Y = Y(m:-1:1,:);
% 
% Uplot = reshape(Uo,n,m);
% Uplot = Uplot';
% Uplot = Uplot(m:-1:1,:);
% Vplot = reshape(Vo,n,m);
% Vplot = Vplot';
% Vplot = Vplot(m:-1:1,:);
% 
% Nx = 5;
% Ny = 4;
% ep = (0.5 + 0.5*JacobiGL(0,0,Nx))';
% et = (0.5 + 0.5*JacobiGL(0,0,Ny));
% 
% x0dg = 1.0;
% y0dg = 0.65;
% lxdg = 0.3;
% lydg = 0.3;
% Jx = 1/lxdg;
% Jy = 1/lydg;
% 
% XDGl = kron(ones(Ny+1,1),ep);
% YDGl = kron(et,ones(1,Nx+1));
% 
% XDGg = lxdg*XDGl + x0dg;
% YDGg = lydg*YDGl + y0dg;
% XDGgt0 = XDGg;
% YDGgt0 = YDGg;
% 
% % Derivative matrices
% XDD = dmatrix(ep,Nx);
% YDD = dmatrix(et,Ny);
% 
% %3 Velocity field interpolated at particle locations
% UDG = interp2(X,Y,Uplot,XDGg,YDGg,'makima');
% VDG = interp2(X,Y,Vplot,XDGg,YDGg,'makima');
% 
% dt=0.1;
% %4. The particles are advected one step in time
% XDGgt1 = XDGg + dt*UDG;
% YDGgt1 = YDGg + dt*VDG;
% NT=1;   
% 
% UDG = interp2(X,Y,Uplot,XDGgt1,YDGg,'makima');
% VDG = interp2(X,Y,Vplot,XDGgt1,YDGg,'makima');
% XDGgt2 = XDGgt1 + dt*UDG;
% YDGgt2 = YDGgt1 + dt*VDG;
% 
% %Computing parametrizations
% xb1 = XDGgt1(1,:);   yb1 = YDGgt1(1,:);
% xb2 = XDGgt1(:,end); yb2 = YDGgt1(:,end);
% xb3 = XDGgt1(end,:); yb3 = YDGgt1(end,:);
% xb4 = XDGgt1(:,1);   yb4 = YDGgt1(:,1);
% [ppgx,ppgy,ppgpx,ppgpy] = param(ep,et,xb1,xb2,xb3,xb4,yb1,yb2,yb3,yb4);
% 
% %Mapping to local element
% x1 = [XDGgt1(1,1); YDGgt1(1,1)];
% x2 = [XDGgt1(1,end); YDGgt1(1,end)];
% x3 = [XDGgt1(end,end); YDGgt1(end,end)];
% x4 = [XDGgt1(end,1); YDGgt1(end,1)];
% [XDGlt1,YDGlt1] = map_to_local(ep,et,XDGgt1,YDGgt1,ppgx,ppgy,ppgpx,ppgpy,x1,x2,x3,x4);
% 
% %Interpolating
% XDGgt2q = lagrange_2D_interpoaltion(ep,et,XDGlt1,YDGlt1,XDGgt2,Nx,Ny);
% YDGgt2q = lagrange_2D_interpoaltion(ep,et,XDGlt1,YDGlt1,YDGgt2,Nx,Ny);
% XDGgt0q = lagrange_2D_interpoaltion(ep,et,XDGlt1,YDGlt1,XDGgt0,Nx,Ny);
% YDGgt0q = lagrange_2D_interpoaltion(ep,et,XDGlt1,YDGlt1,YDGgt0,Nx,Ny);
% 
% %Computing Green-Cauchy tensor and FTLE
% FTLEf = zeros(size(XDGgt0));
% FTLEb = zeros(size(XDGgt0));
% 
% [dxdep,dydep,dxdet,dydet] = derv_epet(Nx,Ny,ep,et,ppgx,ppgy,ppgpx,ppgpy,x1,x2,x3,x4);
%     
% dXDGgt2dep = (XDD*XDGgt2')';
% dYDGgt2dep = (XDD*YDGgt2')';
% dXDGgt2det = YDD*XDGgt2;
% dYDGgt2det = YDD*YDGgt2;
% 
% dphixdx = (1/Jx)*(dXDGgt2dep.*dydet - dXDGgt2det.*dydep);
% dphixdy = (1/Jy)*(dXDGgt2det.*dxdep - dXDGgt2dep.*dxdet);
% dphiydx = (1/Jx)*(dYDGgt2dep.*dydet - dYDGgt2det.*dydep);
% dphiydy = (1/Jy)*(dYDGgt2det.*dxdep - dYDGgt2dep.*dxdet);
% 
% for j=1:Nx+1
%     for i=1:Ny+1
%         gcT = [dphixdx(i,j) dphiydx(i,j);dphixdy(i,j) dphiydy(i,j)];
%         gcT = gcT'*gcT;
%         FTLEf(i,j) = (1/(dt*NT))*log(sqrt(eigs(gcT,1)));
%     end
% end
% 
% dXDGgt2dep = (XDD*XDGgt0')';
% dYDGgt2dep = (XDD*YDGgt0')';
% dXDGgt2det = YDD*XDGgt0;
% dYDGgt2det = YDD*YDGgt0;
% 
% dphixdx = (1/Jx)*(dXDGgt2dep.*dydet - dXDGgt2det.*dydep);
% dphixdy = (1/Jy)*(dXDGgt2det.*dxdep - dXDGgt2dep.*dxdet);
% dphiydx = (1/Jx)*(dYDGgt2dep.*dydet - dYDGgt2det.*dydep);
% dphiydy = (1/Jy)*(dYDGgt2det.*dxdep - dYDGgt2dep.*dxdet);
% 
% for j=1:Nx+1
%     for i=1:Ny+1
%         gcT = [dphixdx(i,j) dphiydx(i,j);dphixdy(i,j) dphiydy(i,j)];
%         gcT = gcT'*gcT;
%         FTLEb(i,j) = (1/(dt*NT))*log(sqrt(eigs(gcT,1)));
%     end
% end
% 
% 
% 
% % figure(2)
% % hold on
% % plot3(XDGgt0,YDGgt0,10*ones(size(XDGg)),'*','Color','b')
% % plot3(XDGgt1,YDGgt1,10*ones(size(XDGg)),'*','Color','k')
% % plot3(XDGgt2,YDGgt2,10*ones(size(XDGg)),'*','Color','g')
% % view(0,90);
% 
% % figure(2)
% % hold on
% % plot3(XDGl,YDGl,10*ones(size(XDGg)),'*','Color','b')
% % plot3(XDGlt1,YDGlt1,10*ones(size(XDGg)),'*','Color','k')
% % view(0,90);
% 
% figure(2)
% hold on
% plot3(XDGgt2,YDGgt2,10*ones(size(XDGg)),'*','Color','b')
% plot3(XDGgt2q,YDGgt2q,10*ones(size(XDGg)),'o','Color','k')
% view(0,90);
% 
% figure(3)
% hold on
% plot3(XDGgt0,YDGgt0,10*ones(size(XDGg)),'*','Color','b')
% plot3(XDGgt0q,YDGgt0q,10*ones(size(XDGg)),'o','Color','k')
% view(0,90);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Initial_condition_60_251_91_Re200.mat

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

Uplot = reshape(Uo,n,m);
Uplot = Uplot';
Uplot = Uplot(m:-1:1,:);
Vplot = reshape(Vo,n,m);
Vplot = Vplot';
Vplot = Vplot(m:-1:1,:);

Nx = 30;
Ny = 30;
ep = (0.5 + 0.5*JacobiGL(0,0,Nx))';
et = (0.5 + 0.5*JacobiGL(0,0,Ny));

x0dg = 0.8;
y0dg = 0.0;
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
% FTLE
[FTLEf,FTLEb,XDGgt2,YDGgt2,XDGgt1,YDGgt1,XDGgt0,YDGgt0] = ...
    FTLE_dgelem(X,Y,U,V,ep,et,XDGgt0,YDGgt0,dt,NT,Nx,Ny,Jx,Jy,XDD,YDD);


figure(2)

% Plotting U velocity
ax1 = subplot(2,2,1);
surf(X,Y,U);
view(0,90);
% caxis([-0.4 0.8])
colormap(ax1,jet)
shading interp
axis equal
axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -10 10]);
colorbar
title('U'); xlabel('X'); ylabel('Y');
drawnow
hold off
    
% Plotting V velocity
ax2 = subplot(2,2,2);
surf(X,Y,V);
view(0,90);
% caxis([-0.4 0.4])
colormap(ax2,jet)
shading interp
axis equal
axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -10 10]);
colorbar
title('V'); xlabel('X'); ylabel('Y');
drawnow

% Plotting FTLE-f
ax3 = subplot(2,2,3);
surf(XDGgt1,YDGgt1,FTLEf);
view(0,90);
% caxis([-0.4 0.4])
shading interp
axis equal
axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -1000 1000]);
colorbar
title('FTLE-f'); xlabel('X'); ylabel('Y');
drawnow

% Plotting FTLE-b
ax4 = subplot(2,2,4);
surf(XDGgt1,YDGgt1,FTLEb);
view(0,90);
% caxis([-0.4 0.4])
shading interp
axis equal
axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -1000 1000]);
colorbar
title('FTLE-b'); xlabel('X'); ylabel('Y');
drawnow


