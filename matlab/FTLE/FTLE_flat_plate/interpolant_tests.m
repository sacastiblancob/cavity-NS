
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
% Nx = 4;
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
% %3 Velocity field interpolated at particle locations
% UDG = interp2(X,Y,Uplot,XDGg,YDGg,'makima');
% VDG = interp2(X,Y,Vplot,XDGg,YDGg,'makima');
% 
% dt=0.1;
% %4. The particles are advected one step in time
% XDGgt1 = XDGg + dt*UDG;
% YDGgt1 = YDGg + dt*VDG;
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
% %INTERPOLATION
% [llai,lmbj] = cross_lagrange(ep,XDGlt1,et,YDGlt1);
% Ipqinv = kron(lmbj,llai);
% psiqx = reshape(XDGgt2',1,[])';
% psipx = Ipqinv*psiqx;
% XDGgt2p = reshape(psipx,Nx+1,[])';
% 
% psiqy = reshape(YDGgt2',1,[])';
% psipy = Ipqinv*psiqy;
% YDGgt2p = reshape(psipy,Nx+1,[])';
% 
% [phipx1,phipy1] = interpol(ep,et,XDGlt1,YDGlt1,XDGgt2,YDGgt2);
% 
% [XDGl2,YDGl2] = map_to_local(ep,et,phipx1,phipy1,ppgx,ppgy,ppgpx,ppgpy,x1,x2,x3,x4);
% 
% % figure(2)
% % hold on
% % plot3(XDGgt2,YDGgt2,10*ones(size(XDGg)),'*','Color','b')
% % % plot3(XDGgt2p,YDGgt2p,10*ones(size(XDGg)),'*','Color','k')
% % plot3(phipx1,phipy1,10*ones(size(XDGg)),'*','Color','g')
% % view(0,90);
% 
% figure(2)
% hold on
% plot3(XDGl,YDGl,10*ones(size(XDGg)),'*','Color','b')
% plot3(XDGl2,YDGl2,10*ones(size(XDGg)),'*','Color','g')
% view(0,90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % %Interpolation in one dimension
% % Nx = 10;
% % ep = (0.5 + 0.5*JacobiGL(0,0,Nx))';
% % yep = ep.*ep;
% % ep1 = 0.3;
% % inep = zeros(size(ep));
% % for i=1:Nx+1
% %     inep(i) = yep(i)*lagrange_point(ep1,ep,i);
% % end
% % yinep = sum(inep);
% 
% %Interpolation in two dimensions
% Nx = 5;
% Ny = 5;
% ep = (0.5 + 0.5*JacobiGL(0,0,Nx))';
% et = (0.5 + 0.5*JacobiGL(0,0,Ny));
% yep = ep.*ep;
% 
% X = kron(ones(Ny+1,1),ep);
% Y = kron(et,ones(1,Nx+1));
% 
% F = X.*X + X.*Y + Y.*Y;
% 
% X2 = X;
% Y2 = Y;
% for j=2:Nx
%     for i=2:Ny
%         X2(i,j) = X2(i,j) + 0.001*i*i;
%         Y2(i,j) = Y2(i,j) + 0.001*j*i;
% %         X2(i,j) = X2(i,j)+0.0001;
% %         Y2(i,j) = Y2(i,j)+0.0001;
%     end
% end
% F1 = X2.*X2 + X2.*Y2 + Y2.*Y2;
% 
% ep1 = ep(2);
% et1 = et(3);
% FR = ep1.*ep1 + ep1.*et1 + et1.*et1;
% 
% XDGl = kron(ones(Ny+1,1),ep);
% YDGl = kron(et,ones(1,Nx+1));
% 
% %INTERPOLATION
% 
% % %Naive way
% % f = zeros(Ny+1,Nx+1);
% % for j=1:Nx+1
% %     for i=1:Ny+1
% %         [i j lagrange_point(ep1,XDGl(i,:),j)]
% %         f(i,j) = F(i,j)*lagrange_point(ep1,XDGl(i,:),j)*lagrange_point(et1,YDGl(:,j),i);
% %     end
% % end
% % F2 = sum(sum(f))
% 
% %Naive way
% f = zeros(Ny+1,Nx+1);
% for j=1:Nx+1
%     for i=1:Ny+1
% %         [i j lagrange_point(ep1,X2(i,:),j)]
%         f(i,j) = F1(i,j)*lagrange_point(ep1,X2(i,:),j)*lagrange_point(et1,Y2(:,j),i);
%     end
% end
% F2 = sum(sum(f));
% 
% % Vectorial way
% [llai,lmbj] = cross_lagrange(ep1,X2,et1,Y2,Nx,Ny);
% psip = reshape(F1,1,[]);
% ll = reshape(llai.*lmbj,1,[])';
% fi = psip*ll;
% 
% 
% % I = zeros((Nx+1)*(Ny+1));
% % k=1;
% % for j=1:(Nx+1)
% %     llai = cross_lagrange_ep(ep(j),X2,Nx,Ny);
% %     for i=1:(Ny+1)
% %         lmbj = cross_lagrange_et(et(i),Y2,Nx,Ny);
% %         I(:,k) = reshape(llai.*lmbj,1,[])';
% %         k=k+1;
% %     end
% % end
% 
% ILA = zeros(Ny+1,(Nx+1)*(Ny+1));
% ILB = zeros((Nx+1)*(Ny+1),Nx+1);
% k=1;
% for j=1:(Nx+1)
%     ILA(:,(j-1)*(Nx+1)+1:j*(Nx+1)) = cross_lagrange_ep(ep(j),X2,Nx,Ny);
% end
% for i=1:(Ny+1)
%     ILB((i-1)*(Ny+1)+1:i*(Ny+1),:) = cross_lagrange_et(et(i),Y2,Nx,Ny);
% end
% 
% psip = reshape(F1,1,[]);
% FI = zeros(Ny+1,Nx+1);
% for j=1:Nx+1
%     for i=1:Ny+1
%         ll = reshape(ILA(:,(j-1)*(Nx+1)+1:j*(Nx+1)).*ILB((i-1)*(Ny+1)+1:i*(Ny+1),:),1,[])';
% %         fi = psip*ll;
%         FI(i,j) = psip*ll;
%     end
% end
%     
% 
% 
% 
% % [liam,ljbl] = cross_lagrange2(ep,flip(X2,2),et,flip(Y2,1),Nx,Ny);
% % I = kron(liam,ljbl);
% % psip = reshape(F1',1,[])';
% % psi = I*psip;
% % FI = reshape(psi,Ny+1,[])
% 
% 
% 
% figure(2)
% hold on
% plot(X,Y,'*','Color','b')
% plot(X2,Y2,'*','Color','k')
% view(0,90);
% 
% % FI = lagrange_2D_interpoaltion(ep,et,X2,Y2,F1,Nx,Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Improving Interpolation and things

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

Nx = 12;
Ny = 12;
ep = (0.5 + 0.5*JacobiGL(0,0,Nx))';
et = (0.5 + 0.5*JacobiGL(0,0,Ny));

% Derivative matrices
XDD = dmatrix(ep,Nx);
YDD = dmatrix(et,Ny);

x0dg = 1.0;
y0dg = 0.65;
lxdg = 0.3;
lydg = 0.3;
Jx = 1/lxdg;
Jy = 1/lydg;

XDGl = kron(ones(Ny+1,1),ep);
YDGl = kron(et,ones(1,Nx+1));

XDGg = lxdg*XDGl + x0dg;
YDGg = lydg*YDGl + y0dg;
XDGgt0 = XDGg;
YDGgt0 = YDGg;

%3 Velocity field interpolated at particle locations
UDG = interp2(X,Y,Uplot,XDGg,YDGg,'makima');
VDG = interp2(X,Y,Vplot,XDGg,YDGg,'makima');

dt=0.1;
%4. The particles are advected one step in time
XDGgt1 = XDGg + dt*UDG;
YDGgt1 = YDGg + dt*VDG;

%Mapping deformed element to local-master element
[XDGl1,YDGl1] = map_to_local_dg(ep,et,XDGgt1,YDGgt1,XDD,YDD);

%Mapping local-master element to deformed element
[XDGlgt1,YDGlgt1] = map_to_global_dg(ep,et,XDGgt1,YDGgt1);

% figure(1)
% hold on
% plot(XDGg,YDGg,'*','Color','b')
% plot(XDGgt1,YDGgt1,'*','Color','k')
% view(0,90);
% 
% figure(2)
% hold on
% plot(XDGl,YDGl,'o','Color','b')
% plot(XDGl1,YDGl1,'*','Color','k')
% view(0,90);
% 
% figure(3)
% hold on
% plot(XDGgt1,YDGgt1,'*','Color','b')
% plot(XDGlgt1,YDGlgt1,'o','Color','k')
% view(0,90);

% %INTERPOLATION
% F = -XDGl1.*XDGl1.*YDGl1 + XDGl1.*YDGl1 + YDGl1.*YDGl1;
% % F = XDGl1 + YDGl1;
% 
% X2=XDGl1;
% Y2=YDGl1;
% F2=F;
% FI = lagrange_2D_interpolation(ep,et,XDGl1,YDGl1,F,Nx,Ny);
% FR = -XDGl.*XDGl.*YDGl + XDGl.*YDGl + YDGl.*YDGl;
% % FR = XDGl + YDGl;
% 
% % norm(FR-FI,'Inf')
% % 
% % surf(XDGl,YDGl,FI);
% % colormap jet
% % hold on
% % s2 = surf(XDGl,YDGl,FR+1);

% %%% Test of speed
% 
% tic
% for i=1:1000
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
% end
% toc
% % Elapsed time is 17.543590 seconds -> 1000 iterations.
% 
% tic
% for i=1:1000
% [XDGl1,YDGl1] = map_to_local_dg(ep,et,XDGgt1,YDGgt1,XDD,YDD);
% end
% toc
% % Elapsed time is 2.802694 seconds -> 1000 iterations.
% %Way better with lagrange interpolant
% 
% 
% % figure(2)
% % hold on
% % plot(XDGl,YDGl,'o','Color','b')
% % plot(XDGl1,YDGl1,'*','Color','k')
% % plot(XDGlt1,YDGlt1,'*','Color','g')
% % view(0,90);

tic
for i=1:1000
%NAIVE WAY
%Computing parametrizations
xb1 = XDGgt1(1,:);   yb1 = YDGgt1(1,:);
xb2 = XDGgt1(:,end); yb2 = YDGgt1(:,end);
xb3 = XDGgt1(end,:); yb3 = YDGgt1(end,:);
xb4 = XDGgt1(:,1);   yb4 = YDGgt1(:,1);
[ppgx,ppgy,ppgpx,ppgpy] = param(ep,et,xb1,xb2,xb3,xb4,yb1,yb2,yb3,yb4);
x1 = [XDGgt1(1,1); YDGgt1(1,1)];
x2 = [XDGgt1(1,end); YDGgt1(1,end)];
x3 = [XDGgt1(end,end); YDGgt1(end,end)];
x4 = [XDGgt1(end,1); YDGgt1(end,1)];

[dxdep,dydep,dxdet,dydet] = derv_epet(Nx,Ny,ep,et,ppgx,ppgy,ppgpx,ppgpy,x1,x2,x3,x4);
end
toc
tic
for i=1:1000
%BEST WAY
[dxdep1,dydep1,dxdet1,dydet1] = derv_epet_dg(ep,et,XDGgt1,YDGgt1,XDD,YDD);
end
toc




