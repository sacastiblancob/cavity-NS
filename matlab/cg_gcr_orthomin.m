%
% This code is based on algorithms for Conjugate Gradient (CG) Conjugate
% Residual (CR) Generalized Conjugate Residual (GCR) ORTHODIR, and others
% of the same kind.
%
% These are complemented with One-dimensional Projection Methods in order
% to stablish relationships.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Matrix A, b, and exact solution
A = [2 1;1 1];
% A = [2 -5;1 1];
b = [0;-1];
xe = A\b;  %exact solution
ni = 10;    %Number of iterations

% h = [1 3 1;1 3 1;2 6 2;2 6 2;3 8 3;3 8 3];
% A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
% A(5,2) = -10;
% A(1,5) = -4;
% n = length(A);
% b = ones(n,1);
% xe = A\b;  %exact solution
% [Av,Ar,Ac] = full2csc(A);
% m=3;
% ni = 10;    %Number of iterations
% tol = 1E-10;
% niter = 50;

%Initial guess
x = [0;-1];
% x = b;

% %GMRES
% [xgm,tgm,res] = csc_gmres(Av,Ar,Ac,b,x,m,niter,tol);

%First residuals and A*r
r = b - A*x;
p = A*r;

rs = r; ps = p; xs = x;        %steepest descend (s)
rm = r; pm = p; xm = x;        %Minimal Residual (m)
rr = r; pr = p; xr = x;        %Residual Norm Steepest Descend (r)
rc = r; pc = r; xc = x;        %Conjugate Gradient (c)
rn = r; pn = r; xn = x;        %Conjugate Residual (n)
rg = r; pg = p; xg = x;        %Generalized Conjugate Residual (g)
ro = r; po = r; xo = x;        %Orthodir (o)

%For storing process
%steepest descend
Rs = zeros(length(b),ni+1); Ps = zeros(length(b),ni+1);
Xs = zeros(length(b),ni+1); Es = zeros(length(b),ni+1); As = zeros(1,ni+1);
Rs(:,1) = rs; Ps(:,1) = ps; Xs(:,1) = xs; Es(:,1) = xe - xs;

%Minimal Residual
Rm = zeros(length(b),ni+1); Pm = zeros(length(b),ni+1);
Xm = zeros(length(b),ni+1); Em = zeros(length(b),ni+1); Am = zeros(1,ni+1);
Rm(:,1) = rm; Pm(:,1) = pm; Xm(:,1) = xm; Em(:,1) = xe - xm;

%Residual Norm Steepest Descend
Rr = zeros(length(b),ni+1); Vr = zeros(length(b),ni+1);
Wr = zeros(length(b),ni+1); Xr = zeros(length(b),ni+1);
Er = zeros(length(b),ni+1);Ar = zeros(1,ni+1);
Rr(:,1) = rr; Xr(:,1) = xr; Er(:,1) = xe - xr;

%Conjugate Gradient
Rc = zeros(length(b),ni+1); Pc = zeros(length(b),ni+1);
Xc = zeros(length(b),ni+1); Ec = zeros(length(b),ni+1);
Ac = zeros(1,ni+1); Bc = zeros(1,ni+1);
Rc(:,1) = rc; Pc(:,1) = pc; Xc(:,1) = xc; Ec(:,1) = xe - xc;

%Conjugate Residual
Rn = zeros(length(b),ni+1); Pn = zeros(length(b),ni+1);
Xn = zeros(length(b),ni+1); En = zeros(length(b),ni+1);
An = zeros(1,ni+1); Bn = zeros(1,ni+1);
Rn(:,1) = rn; Pn(:,1) = pn; Xn(:,1) = xn; En(:,1) = xe - xn;

%Generalized Conjugate Residual
Rg = zeros(length(b),ni+1); Pg = zeros(length(b),ni+1);
Xg = zeros(length(b),ni+1); Eg = zeros(length(b),ni+1);
Ag = zeros(1,ni+1); Bg = zeros(1,ni+1);
Rg(:,1) = rg; Pg(:,1) = pg; Xg(:,1) = xg; Eg(:,1) = xe - xg;

%Generalized Conjugate Residual
Ro = zeros(length(b),ni+1); Po = zeros(length(b),ni+1);
Xo = zeros(length(b),ni+1); Eo = zeros(length(b),ni+1);
Ao = zeros(1,ni+1); Bo = zeros(1,ni+1);
Ro(:,1) = ro; Po(:,1) = po; Xo(:,1) = xo; Eo(:,1) = xe - xo;


for t=1:ni
    %
    %Steepest Descend
    %
    as = (rs'*rs)/(ps'*rs);
    xs = xs + as*rs;
    rs = rs - as*ps;
    ps = A*rs;
    
    %updating
    Rs(:,t+1) = rs; Ps(:,t+1) = ps; Xs(:,t+1) = xs; Es(:,t+1) = xe - xs;
    As(:,t) = as;

    %
    %Minimal Residual Iteration
    %
    am = (pm'*rm)/(pm'*pm);
    xm = xm + am*rm;
    rm = rm - am*pm;
    pm = A*rm;

    %updating
    Rm(:,t+1) = rm; Pm(:,t+1) = pm; Xm(:,t+1) = xm; Em(:,t+1) = xe - xm;
    Am(:,t) = am;

    %
    %Residual Norm Steepest Descend
    %
    vr = A'*rr;
    wr = A*vr;  ar = (vr'*vr)/(wr'*wr);
    xr = xr + ar*vr;
    rr = rr - ar*wr;
    
    %updating
    Rr(:,t+1) = rr; Vr(:,t+1) = vr; Wr(:,t+1) = wr; Xr(:,t+1) = xr;
    Er(:,t+1) = xe - xr; Ar(:,t) = ar;
    
    %
    %ORTHODIR
    %
    ao = (ro'*(A*po))/((A*po)'*(A*po));
    xo = xo + ao*po;
    ro = ro - ao*(A*po);
    bo = -(((A*A)*ro)'*(A*po))/((A*po)'*(A*po));
    po = ro + bo*po;
    
    %updating
    Ro(:,t+1) = ro; Po(:,t+1) = po; Xo(:,t+1) = xo; Eo(:,t+1) = xe - xo;
    Ao(:,t) = ao; Bo(:,t) = bo;
end

for t=1:length(A)   
    %
    %Conjugate Gradient
    %
    ac = (rc'*rc)/((A*pc)'*pc);
    xc = xc + ac*pc;
    rc0 = rc;
    rc = rc - ac*(A*pc);
    bc = (rc'*rc)/(rc0'*rc0);
    pc = rc + bc*pc;
    
    %updating
    Rc(:,t+1) = rc; Pc(:,t+1) = pc; Xc(:,t+1) = xc; Ec(:,t+1) = xe - xc;
    Ac(:,t) = ac; Bc(:,t) = bc;
    
    %
    %Conjugate Residual
    %
    an = (rn'*(A*rn))/((A*pn)'*(A*pn));
    xn = xn + an*pn;
    rn0 = rn;
    rn = rn - an*(A*pn);
    bn = (rn'*(A*rn))/(rn0'*(A*rn0));
    pn = rn + bn*pn;
    
    %updating
    Rn(:,t+1) = rn; Pn(:,t+1) = pn; Xn(:,t+1) = xn; En(:,t+1) = xe - xn;
    An(:,t) = an; Bn(:,t) = bn;
    
    %
    %Generalized Conjugate Residual
    %
    ag = (rg'*(A*pg))/((A*pg)'*(A*pg));
    xg = xg + ag*pg;
    rg = rg - ag*(A*pg);
    bg = -((A*rg)'*(A*pg))/((A*pg)'*(A*pg));
    pg = rg + bg*pg;
    
    %updating
    Rg(:,t+1) = rg; Pg(:,t+1) = pg; Xg(:,t+1) = xg; Eg(:,t+1) = xe - xg;
    Ag(:,t) = ag; Bg(:,t) = bg;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimized Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xmin = min(min(Xs,Xm),Xr);
Xmax = max(max(Xs,Xm),Xr);

xmin = min(Xmin(1,:))-0.2;
xmax = max(Xmax(1,:))+0.2;
if (xmax - xmin)>(4*norm(xe))
    xmin = -3*norm(xe);
    xmax = 3*norm(xe);
end
ymin = min(Xmin(2,:))-0.2;
ymax = max(Xmax(2,:))+0.2;
if (ymax - ymin)>(4*norm(xe))
    ymin = -3*norm(xe);
    ymax = 3*norm(xe);
end
[X,Y] = meshgrid(xmin-0.1:0.1:xmax+0.1,ymin-0.1:0.1:ymax+0.1);
siz = size(X);
Fs = zeros(siz);
Fm = zeros(siz);
F = zeros(siz);
% F2 = zeros(siz);

%F which minimizes Steepest Descend
for i=1:siz(1)
    for j=1:siz(2)
        Fs(i,j) = (A*([X(i,j);Y(i,j)]-xe))'*([X(i,j);Y(i,j)]-xe);
        Fm(i,j) = (b - A*([X(i,j);Y(i,j)]))'*(b - A*([X(i,j);Y(i,j)]));
        F(i,j) = 0.5*([X(i,j) Y(i,j)]*(A*([X(i,j);Y(i,j)]))) - b'*[X(i,j);Y(i,j)];
    end
end
% F2 = Fm - norm(b);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [0 0 20 25]);
% figure(2)
% % plot of X0
% plot([ xe(1,1) xe(1,1)], [xe(2,1) xe(2,1)],'--','Color',[0.7 0.7 0.7]);
% hold all
% plot([ xe(1,1) xe(1,1)], [xe(2,1) xe(2,1)],'-.','Color',[0.7 0.7 0.7]);
% plot([ Xs(1,1) Xs(1,1)], [Xs(2,1) Xs(2,1)],'b');
% % figure(3)
% plot([ Xm(1,1) Xm(1,1)], [Xm(2,1) Xm(2,1)],'r');
% % figure(4)
% plot([ Xr(1,1) Xr(1,1)], [Xr(2,1) Xr(2,1)],'g');
% plot([ Xc(1,1) Xc(1,1)], [Xc(2,1) Xc(2,1)],'--','Color','b')
% plot([ Xn(1,1) Xn(1,1)], [Xn(2,1) Xn(2,1)],'--','Color','r');
% plot([ Xg(1,1) Xg(1,1)], [Xg(2,1) Xg(2,1)],'--','Color','g');
% plot([ Xo(1,1) Xo(1,1)], [Xo(2,1) Xo(2,1)],'--','Color','k');
% 
% figure(2)
% legend({'F-SPD','F-NSPD','SD','MRI','RNSD','CG','CR','GCR','OR'},'AutoUpdate','off','Location','northwest')
% grid on
% axis([xmin-0.2 xmax+0.2 ymin-0.2 ymax+0.2])
% title('All 1D Projection Methods & Conjugate Methods')
% 
% %plot of Minimized Functions
% figure(2)
% cFs = [0.7 0.7 0.7];
% contour(X,Y,Fs,'--','Levels',max(max(Fs))/25)
% % contour(X,Y,F,'--','Levels',max(max(F))/50)
% % contour(X,Y,F2,'*','Levels',max(max(F2))/50)
% colormap(cFs)
% axis equal
% % figure(3)
% % cFm = [0.3 0.3 0.3];
% contour(X,Y,Fm,'-.','Levels',max(max(Fm))/25)
% % contour(X,Y,F2,'*','Levels',max(max(F2))/50)
% % colormap(cFm)
% axis equal
% % figure(4)
% % cFm = [0.3 0.3 0.3];
% % contour(X,Y,Fm,'-.')
% % colormap(cFm)
% % axis equal
% 
% 
% % annotation('textbox',[.9 .5 .4 .2], ...
% %     'String','gray -- , F of SD, gray -., F of MRI and RNSD','EdgeColor','none')
% 
% % figure(3)
% % legend({'MRI','RNSD'},'AutoUpdate','off','Location','northwest')
% % grid on
% % axis([xmin-0.2 xmax+0.2 ymin-0.2 ymax+0.2])
% % title('Min F - Minimal Residual Iteration - Residual Norm SD')
% 
% % figure(4)
% % % legend({'SD'},'AutoUpdate','off','Location','northwest')
% % grid on
% % Xmin = min(min(Xs,Xm),Xr);
% % Xmax = max(max(Xs,Xm),Xr);
% % axis([min(Xmin(1,:))-0.2 max(Xmax(1,:))+0.2 min(Xmin(2,:))-0.2 max(Xmax(2,:))+0.2])
% % title('Residual Norm Steepest Descend')
% 
% ARs = As.*Rs; ARm = Am.*Rm; ARr = Ar.*Rr; ARC = Ac.*Rc; ARn = An.*Rn;
% ARg = Ag.*Rg; ARo = Ao.*Ro;
% 
% for t=1:ni
%     
%     %Steepest Descend
%     figure(2)
%     plot([ Xs(1,t) Xs(1,t+1)], [Xs(2,t) Xs(2,t+1)],'b');
%     
%     %Minimal Residual Iteration
% %     figure(3)
%     plot([ Xm(1,t) Xm(1,t+1)], [Xm(2,t) Xm(2,t+1)],'r');
% 
% %     figure(4)
%     %Residual Norm Steepest Descend
%     plot([ Xr(1,t) Xr(1,t+1)], [Xr(2,t) Xr(2,t+1)],'g');
%     
%     %ORTHODIR
%     plot([ Xo(1,t) Xo(1,t+1)], [Xo(2,t) Xo(2,t+1)],'--','Color','k');
% 
% end
% 
% for t=1:length(A)
%     
%     %Conjugate Gradient
%     figure(2)
%     plot([ Xc(1,t) Xc(1,t+1)], [Xc(2,t) Xc(2,t+1)],'--','Color','b');
%     
%     %Conjugate Residual
% %     figure(3)
%     plot([ Xn(1,t) Xn(1,t+1)], [Xn(2,t) Xn(2,t+1)],'--','Color','r');
% 
% %     figure(4)
%     %Generalized Conjugate Residual
%     plot([ Xg(1,t) Xg(1,t+1)], [Xg(2,t) Xg(2,t+1)],'--','Color','g');
%     
% end
% 
% figure(2)
% axis([xmin-0.2 xmax+0.2 ymin-0.2 ymax+0.2])
% 
% % figure(3)
% % axis([xmin-0.2 xmax+0.2 ymin-0.2 ymax+0.2])
% 
% % figure(5)
% % subplot(3,1,1)
% % surf(X,Y,Fs)
% % subplot(3,1,2)
% % surf(X,Y,F2)
% % subplot(3,1,3)
% % surf(X,Y,F)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of SD-CG and RNSD-GCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xmin = min(Xs,Xr);
Xmax = max(Xs,Xr);

xmin = min(Xmin(1,:))-0.2;
xmax = max(Xmax(1,:))+0.2;
if (xmax - xmin)>(4*norm(xe))
    xmin = -2*norm(xe);
    xmax = 2*norm(xe);
end
ymin = min(Xmin(2,:))-0.2;
ymax = max(Xmax(2,:))+0.2;
if (ymax - ymin)>(4*norm(xe))
    ymin = -2*norm(xe);
    ymax = 2*norm(xe);
end

figure(4)
subplot(1,2,1)
% plot of X0
plot([ xe(1,1) xe(1,1)], [xe(2,1) xe(2,1)],'--','Color',[0.7 0.7 0.7]);
hold all
plot([ Xs(1,1) Xs(1,1)], [Xs(2,1) Xs(2,1)],'b');
plot([ Xc(1,1) Xc(1,1)], [Xc(2,1) Xc(2,1)],'--','Color','b')

legend({'F','SD','CG'},'AutoUpdate','off','Location','northwest')
grid on
axis([xmin-0.2 xmax+0.2 ymin-0.2 ymax+0.2])
title('Steepest Descent (SD) vs Conjugate Gradient (CG)')

cFs = [0.7 0.7 0.7];
contour(X,Y,Fs,'--','Levels',max(max(Fs))/25)
colormap(cFs)

for t=1:ni
    
    %Steepest Descend
    plot([ Xs(1,t) Xs(1,t+1)], [Xs(2,t) Xs(2,t+1)],'b');
    
end

for t=1:length(A)
    
    %Conjugate Gradient
    plot([ Xc(1,t) Xc(1,t+1)], [Xc(2,t) Xc(2,t+1)],'--','Color','b');
    
end

axis equal
axis([xmin xmax ymin ymax])

%-------------------------------------------------------------------------%

figure(4)
subplot(1,2,2)
% plot of X0
plot([ xe(1,1) xe(1,1)], [xe(2,1) xe(2,1)],'--','Color',[0.7 0.7 0.7]);
hold all
plot([ Xr(1,1) Xr(1,1)], [Xr(2,1) Xr(2,1)],'r');
plot([ Xg(1,1) Xg(1,1)], [Xg(2,1) Xg(2,1)],'--','Color','r');

legend({'F','RNSD','GCR'},'AutoUpdate','off','Location','northwest')
grid on
axis([xmin-0.2 xmax+0.2 ymin-0.2 ymax+0.2])
title('Residual Norm Steepest Descend (RNSD) vs Generalized Conjugate Residual (GCR)')

cFm = [0.7 0.7 0.7];
contour(X,Y,Fm,'--','Levels',max(max(Fs))/25)
colormap(cFm)

for t=1:ni
    
    %Residual Norm Steepest Descend
    plot([ Xr(1,t) Xr(1,t+1)], [Xr(2,t) Xr(2,t+1)],'r');
    
end

for t=1:length(A)
    
    %Generalized Conjugate Residual
    plot([ Xg(1,t) Xg(1,t+1)], [Xg(2,t) Xg(2,t+1)],'--','Color','r');
    
end

axis equal
axis([xmin xmax ymin ymax])







