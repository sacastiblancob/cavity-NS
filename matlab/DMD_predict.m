%% Loading Data
% load field_results_281_79_Re300.mat
load field_results_41-101_281_79_Re300.mat
mus = size(MU);
tt = mus(2);
tprint = 10;    %how often time-steps information was taken
umax = max(max(MU)); umin = min(min(MU)); vmax = max(max(MV));
vmin = min(min(MV)); pmax = max(max(MP)); pmin = min(min(MP));

%%
% % Plotting data
% bstep = 5;
% for t=1:bstep:tt+1
% 
%     Uplot = reshape(MU(:,t),n,m);
%     Uplot = Uplot';
%     Uplot = Uplot(m:-1:1,:);
%     Vplot = reshape(MV(:,t),n,m);
%     Vplot = Vplot';
%     Vplot = Vplot(m:-1:1,:);
%     Pplot = reshape(MP(:,t),n,m);
%     Pplot = Pplot';
%     Pplot = Pplot(m:-1:1,:);
% 
% %     Plotting U velocity
%     ax1 = subplot(3,1,1);
%     surf(X,Y,Uplot); 
%     axis([xmin xmax ymin ymax umin-0.1 umax+0.1]);
%     view(0,90);
% % %     caxis([-0.4 1])
%     colormap(ax1,jet)
%     shading interp
%     axis equal
%     colorbar
%     title(strcat('U, Flat Plate of ',string(lp),', Re=',string(round(Re)),', Time =',string(round((t-1)*dt*10,1)))); xlabel('X'); ylabel('Y');
%     drawnow
%         
% %     Plotting V velocity
%     ax2 = subplot(3,1,2);
%     surf(X,Y,Vplot); 
%     axis([xmin xmax ymin ymax vmin-0.1 vmax+0.1]);
%     view(0,90);
% % %     caxis([-0.4 0.1])
%     colormap(ax2,jet)
%     shading interp
%     axis equal
%     colorbar
%     title('V'); xlabel('X'); ylabel('Y');
%     drawnow
%     
% %     Plotting pressure
%     ax3 = subplot(3,1,3);
%     surf(X,Y,Pplot); 
%     axis([xmin xmax ymin ymax pmin-10 pmax+10]);
%     view(0,90);
% % %     caxis([-0.4 0.1])
%     colormap(ax3,gray)
%     shading interp
%     axis equal
%     colorbar
%     title('P'); xlabel('X'); ylabel('Y');
%     drawnow
% 
% end

%% Splitting Data 1/3 - 2/3
pp = 2/3;
tc = round(tt*pp);
tv = tt - tc;
tsamp = 1;
DAT = MV;
DATcal = DAT(:,1:tsamp:tc);
DATval = DAT(:,tc+1:tsamp:tt);

%% CALIBRATION

% COMPUTING SVD
sizd = size(DATcal);
j = sizd(2);
XX = DATcal(:,1:j-1);
XX2 = DATcal(:,2:j);
[U,S,V] = svd(XX,'econ');
ds = diag(S);
% semilogy(ds(1:50))
% semilogy(ds)

% COMPUTING DMD
orcut = 0;
r = sum(ds>=orcut);
Uu = U(:,1:r);
Ss = S(1:r,1:r);
Vv = V(:,1:r);
Si = diag(1./diag(Ss));
At = Uu'*XX2*Vv*Si;
[W, eigs] = eig(At);
Phi = XX2*Vv*Si*W;

%% Prediction
MU0 = DATcal(:,end);
b0 = Phi'*MU0;

DATpre = zeros(mus(1),tv);
DT = tprint*tsamp*dt;
eigs = diag(diag(eigs));
for ti = 1:tv
    DATpre(:,ti) = Phi*(eigs.^ti)*b0;
end
DATpre = real(DATpre);

%% Validation
normerror = zeros(1,tv);
norminf = zeros(1,tv);
for i=1:tv
    normerror(i) = norm(DATpre(:,i) - DATval(:,i));
    norminf(i) = norm(DATpre(:,i) - DATval(:,i),'Inf');
end
figure(1)
subplot(1,2,1)
semilogy(normerror)
subplot(1,2,2)
semilogy(norminf)

%% Plotting solutions
bstep = 5;
zmax1 = max(max(DATpre)); zmin1 = min(min(DATpre));
zmax2 = max(max(DATval)); zmin2 = min(min(DATval));
figure(2)
for t=1:bstep:tv
% for t=1

    Uplot = reshape(DATpre(:,t),n,m);
    Uplot = Uplot';
    Uplot = Uplot(m:-1:1,:);
    Vplot = reshape(DATval(:,t),n,m);
    Vplot = Vplot';
    Vplot = Vplot(m:-1:1,:);

%     Plotting U velocity
    ax1 = subplot(2,1,1);
    surf(X,Y,Uplot); 
    axis([xmin xmax ymin ymax zmin1-0.1 zmax1+0.1]);
    view(0,90);
    caxis([zmin2-0.1 zmax2+0.1])
    colormap(ax1,jet)
    shading interp
    axis equal
    colorbar
    title('Prediction'); xlabel('X'); ylabel('Y');
    drawnow
        
%     Plotting V velocity
    ax2 = subplot(2,1,2);
    surf(X,Y,Vplot); 
    axis([xmin xmax ymin ymax zmin2-0.1 zmax2+0.1]);
    view(0,90);
    caxis([zmin2-0.1 zmax2+0.1])
    colormap(ax2,jet)
    shading interp
    axis equal
    colorbar
    title('Validation'); xlabel('X'); ylabel('Y');
    drawnow
    
end
