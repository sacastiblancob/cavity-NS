%% Loading Data
% load field_results_281_79_Re300.mat
load field_results_FTLE_40-60_251_91_Re200.mat
mus = size(MU);
% tt = mus(2);
tt = 25;
tsamp = 20;    %how often time-steps information was taken
umax = max(max(MU)); umin = min(min(MU)); vmax = max(max(MV));
vmin = min(min(MV)); pmax = max(max(MP)); pmin = min(min(MP));

%%
% % Plotting data
bstep = 1;
for t=1:bstep:tt+1

    Uplot = reshape(MU(:,t),n,m);
    Uplot = Uplot';
    Uplot = Uplot(m:-1:1,:);
    Vplot = reshape(MV(:,t),n,m);
    Vplot = Vplot';
    Vplot = Vplot(m:-1:1,:);
    FTLEF = reshape(MFTLEf(:,t),sizeg(1),sizeg(2));
    FTLEB = reshape(MFTLEb(:,t),sizeg(1),sizeg(2));

%     Plotting U velocity
    ax1 = subplot(2,2,1);
    surf(X,Y,Uplot); 
    axis([xmin xmax ymin ymax umin-0.1 umax+0.1]);
    view(0,90);
    caxis([umin umax])
    colormap(ax1,jet)
    shading interp
    axis equal
    colorbar
    title(strcat('U, Flat Plate of ',string(lp),', Re=',string(round(Re)),', Time =',string(round((t-1)*dt*10,1)))); xlabel('X'); ylabel('Y');
    drawnow
        
%     Plotting V velocity
    ax2 = subplot(2,2,2);
    surf(X,Y,Vplot); 
    axis([xmin xmax ymin ymax vmin-0.1 vmax+0.1]);
    view(0,90);
    caxis([vmin vmax])
    colormap(ax2,jet)
    shading interp
    axis equal
    colorbar
    title('V'); xlabel('X'); ylabel('Y');
    drawnow
    
%     FTLE_F field
    ax3 = subplot(2,2,3);
    surf(XG,YG,FTLEF);
    view(0,90);
    caxis([-1 0.2])
    colormap(ax3,jet)
    shading interp
    axis equal
    axis([xmin xmax ymin ymax -10 10]);
    colorbar
    title('FTLE-f'); xlabel('X'); ylabel('Y');
    drawnow

    %     FTLE_B field
    ax4 = subplot(2,2,4);
    surf(XG,YG,FTLEB);
    view(0,90);
    caxis([-1 0.2])
    colormap(ax4,jet)
    shading interp
    axis equal
    axis([xmin xmax ymin ymax -10 10]);
    colorbar
    title('FTLE-b'); xlabel('X'); ylabel('Y');
    drawnow

end


