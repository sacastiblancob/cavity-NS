load field_results_281_79_Re300.mat

%% Loading Data
DAT = MU;
sizd = size(DAT);
j = sizd(2);
XX = DAT(:,1:j-1);
XX2 = DAT(:,2:j);
[U,S,V] = svd(XX,'econ');
semilogy(diag(S))

%% Plotting data
for t=1:T+1

    Uplot = reshape(MU(:,t),n,m);
    Uplot = Uplot';
    Uplot = Uplot(m:-1:1,:);
    Vplot = reshape(MV(:,t),n,m);
    Vplot = Vplot';
    Vplot = Vplot(m:-1:1,:);
    Pplot = reshape(MP(:,t),n,m);
    Pplot = Pplot';
    Pplot = Pplot(m:-1:1,:);

    % Plotting U velocity
    ax1 = subplot(3,1,1);
    surf(X,Y,Uplot); 
    axis([xmin xmax ymin ymax -1 1.2]);
    view(0,90);
%     caxis([-0.4 1])
    colormap(ax1,jet)
    shading interp
    axis equal
    colorbar
    title(strcat('U, Flat Plate of ',string(lp),', Re=',string(round(Re)),', Time =',string(round((t-1)*dt,1)))); xlabel('X'); ylabel('Y');
    drawnow
        
    % Plotting V velocity
    ax2 = subplot(3,1,2);
    surf(X,Y,Vplot); 
    axis([xmin xmax ymin ymax -1 1]);
    view(0,90);
%     caxis([-0.4 0.1])
    colormap(ax2,jet)
    shading interp
    axis equal
    colorbar
    title('V'); xlabel('X'); ylabel('Y');
    drawnow
    
    % Plotting pressure
    ax3 = subplot(3,1,3);
    surf(X,Y,Pplot); 
    axis([xmin xmax ymin ymax -1000 1000]);
    view(0,90);
%     caxis([-0.4 0.1])
    colormap(ax3,gray)
    shading interp
    axis equal
    colorbar
    title('P'); xlabel('X'); ylabel('Y');
    drawnow

end

%% Compute DMD
r = 21;
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Si = diag(1./diag(S));
At = U'*XX2*V*Si;
[W, eigs] = eig(At);
Phi = XX2*V*Si*W;

%% Plotting DMD

%Number of nodes in each direction
nm = 4;                      %Number of modes to plot

for i=1:nm

% Dplot = reshape(real(Phi(:,i)),n,m);
Dplot = reshape(imag(Phi(:,i)),n,m);
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
title(strcat('U_Modes=',string(i),', Flat Plate of=',string(lp),', P-Re=',string(round(Re)))); xlabel('X'); ylabel('Y');
drawnow

end


