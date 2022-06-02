%% Loading Data
% load field_results_281_79_Re300.mat
load field_results_41-101_281_79_Re300.mat
mus = size(MU);
tt = mus(2);
tsamp = 10;    %how often time-steps information was taken
umax = max(max(MU)); umin = min(min(MU)); vmax = max(max(MV));
vmin = min(min(MV)); pmax = max(max(MP)); pmin = min(min(MP));

% %%
% % % Plotting data
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
%     caxis([umin umax])
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
%     caxis([vmin vmax])
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
%     caxis([pmin pmax])
%     colormap(ax3,gray)
%     shading interp
%     axis equal
%     colorbar
%     title('P'); xlabel('X'); ylabel('Y');
%     drawnow
% 
% end

%% Computing SVD
tsamp = 1;
DAT = MU(:,1:tsamp:tt);
sizd = size(DAT);
j = sizd(2);
XX = DAT(:,1:j-1);
XX2 = DAT(:,2:j);
[U,S,V] = svd(XX,'econ');
ds = diag(S);
% semilogy(ds(1:50))
% semilogy(ds)

%% Compute DMD
orcut = 100;
r = sum(ds>=orcut);
Ss = S(1:r,1:r);
Vv = V(:,1:r);
Uu = U(:,1:r);
Si = diag(1./diag(Ss));
At = Uu'*XX2*Vv*Si;
[W, eigs] = eig(At);
Phi = XX2*Vv*Si*W;

%% Plotting DMD

%Number of nodes in each direction
nm = 4;                      %Number of modes to plot

for i=1:nm

Dplot = reshape(real(Phi(:,i)),n,m);
% Dplot = reshape(imag(Phi(:,i)),n,m);
Dplot = Dplot';

figure(i)
% Plotting U velocity
ax1 = subplot(1,1,1);
surf(X,Y,Dplot); 
axis([xmin xmax ymin ymax -1 1.2]);
view(0,90);
%     caxis([-0.4 1])
colormap(ax1,jet)
shading interp
axis equal
colorbar
title(strcat('Modes=',string(i),', Flat Plate of=',string(lp),', P-Re=',string(round(Re)))); xlabel('X'); ylabel('Y');
drawnow

end


