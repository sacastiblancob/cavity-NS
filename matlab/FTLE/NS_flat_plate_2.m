% clear; clc;
% clear all
% tic
%%
%Solver Burgers 2D, with primitive variables scheme for the nonlinear
%terms.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical paramters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nu = 0.0000392;          %Kinematic viscosity (ISU)
nu = 0.00045;
rho = 1000;         %Density of water, kg/m3 (ISU)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Space parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of nodes in each direction
n = 251;                      %Number of nodes in X direction
m = 91;                      %Number of nodes in Y direction

%Boundaries locations of the domain
xmin = 0.0;       %Left boundary
ymin = -0.4;       %Bottom boundary
xmax = 5.0;       %Top boundary
ymax = 1.4;       %Right boundary

%Coordinates of location of flat plate
xp = 0.8;
yp = 0.5;

%Length of flat plate
lp = 0.3;

%Angle of attack (radians)
alfa = pi/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet values of velocity in X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Left
dirl = 0.3;
% %Rigth
% dirr = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Initial Condition if Any
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save('Initial_condition_60_251_91_Re200.mat','Uo','Vo','P','dt','dirl','n','m','xmin','ymin','xmax','ymax','xp','yp','lp','nu','rho')
% load Initial_condition_20_251_91_Re200.mat
load Initial_condition_40_251_91_Re200.mat
% load Initial_condition_60_251_91_Re200.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Courant Number
CFL = 0.25;                                    %CFL (nonlinear parameter)

% Initial time
to = 0;                                     %Initial time

% Final time
tf = 20;                                  %Final time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing Space Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
pos = 1:n*m;
Mpos = zeros(m,n);
Mpos(m,:) = 1:n;
for i=m-1:-1:1
    Mpos(i,:) = Mpos(i+1,:) + n;
end

xaux = xmin:(xmax-xmin)/10:xmax;
y = [1 0.5 1 0.5 1 0.5 1 0.5 1 0.5 1];
pp = spline(xaux,y);
x2 = (x - xmin);
x2 = x2/max(x2);
x2 = 10*x2;

%Boundary nodes
bound = [1:n,n+1:n:n*m-2*n+1,2*n:n:n*m-n,n*m-n+1:n*m];
bound = sort(bound);
upbound = n*m-n+2:n*m-1;
dobound = 2:n-1;
lebound = n+1:n:n*m-2*n+1;
ribound = 2*n:n:n*m-n;
upboundin = (n*m-2*n+3:n*m-n-2);
doboundin = (n+3:2*n-2);
leboundin = (2*n+2:n:n*m-3*n+2);
riboundin = (3*n-1:n:n*m-2*n-1);
boundint = [Mpos(2,2:n-1), Mpos(m-1,2:n-1), Mpos(3:m-2,2)',Mpos(3:m-2,n-1)'];
boundint = sort(boundint);

%Vecinity of dirichlet boundaries (left and right)
dboundlin = leboundin;
dboundrin = riboundin;
dboundl = lebound;
dboundr = ribound;
dbound = bound;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial condition for velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uo = zeros(m*n,1);
% Uo(lebound) = dirl;
% % Uo(ribound) = dirr;
% Vo = zeros(m*n,1);
% % % Initial P
% P = ones(n*m,1);
% P(lebound) = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Temporal parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dt = (min(dx,dy)*CFL/(max(abs(Uo))));    %Size of the time step
T = 10*round((tf-to)/(10*dt));                      %Number of time steps
t = to+dt:dt:(T*dt)+to;                        %Computational time

%Reynolds Number
% Re = max(Uo)*max([(xmax-xmin) (ymax-ymin)])/nu;
Re = dirl*lp/nu;

disp('Reynolds Number')
fprintf('Re = %.0f \n',Re)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stiffnes linear diffusion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sx and Sy
Sx = ((dt*nu)/(dx^2));  %Linear for x
Sy = ((dt*nu)/(dy^2));  %Linear for y

%constants
ap=(1+2*Sx+2*Sy);
ax=-Sx;
ay=-Sy;

% Computing matrix K for diffusion with Kronecker products
K11 = spdiags(ones(m-2,1),0,(m-2),(m-2));
K12 = spdiags([ax*ones(n-2,1) ap*ones(n-2,1) ax*ones(n-2,1)],[-1 0 1],(n-2),(n-2));
K21 = spdiags(ones(m-2,2),[-1 1],(m-2),(m-2));
K22 = spdiags(ay*ones(n-2,1),0,n-2,n-2);
K1 = kron(K11,K12);
K2 = kron(K21,K22);
K = K1 + K2;
%spy(K)
clear K11 K12 K21 K22 K1 K2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Laplace operator matrix with five points second drivate and Neumann
%boundary conditions (centered difference and ghost nodes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Storing diagonals in one matrix and their positions
% Ensuring 
if(dy>=dx)
    alx = 1;
    aly = (dx^2)/(dy^2);
    alp = -(2 + ((dx^2)/(dy^2))*2);
else
    alx = (dy^2)/(dx^2);
    aly = 1;
    alp =-(2*((dy^2)/(dx^2)) + 2);
end

% %
% % Computing matrix L for Poisson Equation with Kronecker products
% % Standard Matlab

L11 = spdiags(ones(m,1),0,m,m);
h1 = [alx*ones(n,1) alp*ones(n,1) alx*ones(n,1)];
h1(2,3) = alx*2;
h1(n-1,1) = alx*2;
L12 = spdiags(h1,[-1 0 1],n,n);
h2 = ones(m,2);
h2(2,2) = 2;
h2(m-1,1) = 2;
L21 = spdiags(h2,[-1 1],m,m);
L22 = spdiags(aly*ones(n,1),0,n,n);
L1 = kron(L11,L12);
L2 = kron(L21,L22);
L = L1 + L2;
% % %spy(L)
clear L11 L12 L21 L22 L1 L2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Including flate plate modifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%flate plate nodes
%

% % slope of plate
if (abs(alfa)-(pi/2))<=1E-8

    %integer coordinate of corner of table
    npc = round((xp - xmin)/dx);
    mpc = round(((yp-lp/2) - ymin)/dy)+1;
    ipc = (mpc-1)*n + npc;

    %Initializing plate indexes 
    upplate = [];
    doplate = [];

    %number of integers up
    niu = round(lp/dy)+1;
    riplate = zeros(1,niu);
    leplate = zeros(1,niu);
    for k=1:niu
        leplate(k) = ipc + (k-1)*n;
        riplate(k) = leplate(k) + 1;
    end

    for i=leplate
        L(i,:) = 0;
        L(i,i) = alp;
        L(i,i-1) = 2*alx;
        L(i,i+n) = 1*aly;
        L(i,i-n) = 1*aly;

        L(i+1,:) = 0;
        L(i+1,i+1) = alp;
        L(i+1,i+2) = 2*alx;
        L(i+1,i+1+n) = 1*aly;
        L(i+1,i+1-n) = 1*aly;
    end
    
    k=1;
    for i=1:n*m
        if ismember(i,bound)
            continue
        elseif (ismember(i,leplate))
            K(k,:) = 0;
            K(k,k) = ap;
            K(k,k-1) = 2*ax;
            K(k,k+n-2) = 1*ay;
            K(k,k-n+2) = 1*ay;
            k = k + 1;
        elseif (ismember(i,leplate+1))
            K(k,:) = 0;
            K(k,k) = ap;
            K(k,k+1) = 2*ax;
            K(k,k+n-2) = 1*ay;
            K(k,k-n+2) = 1*ay;
            k = k + 1;
        else
            k = k + 1;
        end
    end
    
    bound = sort([bound leplate leplate+1]);
    lebound = sort([lebound leplate+1]);
    ribound = sort([ribound leplate]);
    leboundin = sort([leboundin leplate+2]);
    riboundin = sort([riboundin leplate-1]);

elseif (abs(alfa))<=1E-8

else
    %length X and Y of flat plate
    lxp = lp*sin(alfa);
    lyp = lp*cos(alfa);

    %slope of flat plate
    mp = sin(alfa)/cos(alfa);

%integer
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing minimum singular vector related of L' for regularization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% C. POZRIKIDIS, 2001
% A Note on the Regularization of the Discrete Poisson–Neumann Problem
%

% % % WITH STANDARD MATLAB FUNCTION
%Left and Right hand side Singular Vector matlab
[UL,SL,VL] = svds(L,1,'smallest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL DG MESH

%Number of number of subdomains in X and Y
Nxg = 40;
Nyg = 12;

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
% XG1 = XG; YG1 = YG;
% XG2 = XG; YG2 = YG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOR STORAGE OF RESULTS

j = 1;
rv = 20;

MU = zeros(n*m,round(T/rv)+1);
MV = zeros(n*m,round(T/rv)+1);
MP = zeros(n*m,round(T/rv)+1);
MU(:,1) = Uo;
MV(:,1) = Vo;
MP(:,1) = P;
tt = 2;
hh=2;

sizeg = size(XG);
FTLEF = zeros(size(XG)); FTLEB = zeros(size(XG));
MFTLEf = zeros(sizeg(1)*sizeg(2),round(T/rv)+1);
MFTLEb = zeros(sizeg(1)*sizeg(2),round(T/rv)+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing FTLE fields for the initial condition

%CFL, NT
CFLftle = 3;
maxvel = max(max(abs(Uo),abs(Vo)));
dtftle = CFLftle*(min(dx,dy)/maxvel);
NTf = 40;

%Velocities for interpolation
UI = reshape(Uo,n,m);
UI = UI';
UI = UI(m:-1:1,:);
VI = reshape(Vo,n,m);
VI = VI';
VI = VI(m:-1:1,:);

UI =[UI(:,1) UI zeros(m,1)];
VI =[zeros(m,1) VI zeros(m,1)];
XI =[(xmin-2)*ones(m,1) X (xmax+2)*ones(m,1)];
YI =[Y(:,1) Y Y(:,1)];

for j=1:Nyg
    for i=1:Nxg
        XDGgt0 = XG((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1);
        YDGgt0 = YG((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1);
        
        % FTLE-forward
        Jx = JX(i);
        Jy = JY(j);
        [FTLEf] = ...
            FTLE_dgelem_f(XI,YI,UI,VI,XDGgt0,YDGgt0,dtftle,NTf,Nx,Ny,Jx,Jy,XDD,YDD);

        %Storing results
        FTLEF((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = FTLEf;
        
        % FTLE-backward
        [FTLEb] = ...
            FTLE_dgelem_f(XI,YI,-UI,-VI,XDGgt0,YDGgt0,dtftle,NTf,Nx,Ny,Jx,Jy,XDD,YDD);
        
        %Storing results
        FTLEB((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = FTLEb;

    end
end
MFTLEf(:,1) = reshape(FTLEF,[],1);
MFTLEb(:,1) = reshape(FTLEB,[],1);

%%
%Temporal loop
for time=t
     
    time
    max(Uo)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First fractional step (Non-lineal advection)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Derivatives
    [dudx,dudy,dvdx,dvdy] = diver2D(Uo,Vo,dx,dy,n,m,bound,...
        upboundin,doboundin,leboundin,riboundin);
    
    %Resultado (Explicito)
    Up = Uo - dt*(Uo.*dudx + Vo.*dudy);
    Vp = Vo - dt*(Uo.*dvdx + Vo.*dvdy);
    
%     %Boundary conditions
%     Up(upbound) = 1.0;
%     Up(n*m) = 1.0;
%     Up(n*m-n+1) = 1.0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Second fractional step (Poisson equation)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Derivacion de U y V resultado del primer paso fraccionado en sus respectivas direcciones principales
    [dupdx,dvpdy] = grad2D(Up,Vp,dx,dy,n,m,upbound,dobound,lebound,ribound);
    
    %Vector de mano derecha sin regularizar
    if(dy>=dx)
        rhsp = ((dx^2)*rho/dt)*(dupdx + dvpdy);
    else
        rhsp = ((dy^2)*rho/dt)*(dupdx + dvpdy);
    end
        
    %
    %%% Regularizacion (Posrikidiz, 2001) %%%
    %
    
    % Regularized right hand side
%     rhspr = RM*rhsp;
%     if j<2
        rhtu = rhsp'*UL;
        rhspr = rhsp - rhtu*UL;
%     else
%         rhspr = rhsp;
%     end
    
        
    %%% Neumann Boundary Condition == 0 for Poisson Equation %%%
    %%% They are already included in matrix L structure
    %%% Since Neumann boundary condition == 0, there's nothing to do with
    %%% the right hand side (due to the use of central difference and ghost
    %%% nodes for taking into account Neumann boundary condition)
    
    %
    %%% Solving poisson equation %%%
    %
    
    % % % With Matlab Standard function
    P = L\rhspr;
    j=j+1;
    
    % Computing pressure gradient
    [dpdx,dpdy] = grad2D(P,P,dx,dy,n,m,upbound,dobound,lebound,ribound);
    Upp = Up - (dt/rho)*dpdx;
    Vpp = Vp - (dt/rho)*dpdy;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Third fractional step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Taking internal portion of Upp and Vpp (diffusion RHS) and settling
    % dirichlet boundary conditions in the RHS
    drhsx = zeros((n-2)*(m-2),1);
    drhsy = zeros((n-2)*(m-2),1);
    k=1;
    for i=1:n*m
        if ismember(i,dbound)
            continue
        elseif (ismember(i,dboundlin))
            drhsx(k) = Upp(i) - ax*dirl;
            drhsy(k) = Vpp(i);
            k = k + 1;
        elseif (ismember(i,dboundrin))
%             drhsx(k) = Upp(i) - ax*dirr;
            drhsx(k) = Upp(i);
            drhsy(k) = Vpp(i);
            k = k + 1;
        elseif (ismember(i,leplate))
            drhsx(k) = 0;
            drhsy(k) = 0;
            k = k + 1;
        elseif (ismember(i,leplate+1))
            drhsx(k) = 0;
            drhsy(k) = 0;
            k = k + 1;
        else
            drhsx(k) = Upp(i);
            drhsy(k) = Vpp(i);
            k = k + 1;
        end
    end

    % Solving diffusion step with Standard Matlab solver
    Ud = K\drhsx;
    Vd = K\drhsy;
    

    % Adding boundary conditions
    U = zeros(n*m,1);
    V = zeros(n*m,1);
    k = 1;
    for i=1:n*m
        if ismember(i,dbound)
            if (ismember(i,dboundl))
                U(i) = dirl;
                V(i) = 0;
            elseif (ismember(i,dboundr))
%                 U(i) = dirr;
%                 U(i) = 0.25*(2*U(i-1) + U(i-1-n) + U(i-1+n));
%                 V(i) = 0.25*(2*V(i-1) + V(i-1-n) + V(i-1+n));
                U(i) = U(i-1);
                V(i) = V(i-1);
            end
        else
            U(i) = Ud(k);
            V(i) = Vd(k);
            k = k + 1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Uo and Vo and Storing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Uo = U;
    Vo = V;

    if mod(tt,rv)==0
        MU(:,hh) = Uo;
        MV(:,hh) = Vo;
        MP(:,hh) = P;
%         hh = hh+1;
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Postprocessing and plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Uplot = reshape(Uo,n,m);
    Uplot = Uplot';
    Uplot = Uplot(m:-1:1,:);
    Vplot = reshape(Vo,n,m);
    Vplot = Vplot';
    Vplot = Vplot(m:-1:1,:);
    Pplot = reshape(P-P(n),n,m);
    Pplot = Pplot';
    Pplot = Pplot(m:-1:1,:);

        %Computing LCS - DG
    if mod(tt,rv)==0

%     %CFL, NT
%     CFLftle = 3;
%     maxvel = max(max(abs(Uo),abs(Vo)));
%     dtftle = CFLftle*(min(dx,dy)/maxvel);

    %Velocities for interpolation
    UI =[Uplot(:,1) Uplot zeros(m,1)];
    VI =[zeros(m,1) Vplot zeros(m,1)];

    for j=1:Nyg
        for i=1:Nxg
            XDGgt0 = XG((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1);
            YDGgt0 = YG((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1);
            
            % FTLE-forward
            Jx = JX(i);
            Jy = JY(j);
            [FTLEf] = ...
                FTLE_dgelem_f(XI,YI,UI,VI,XDGgt0,YDGgt0,dtftle,NTf,Nx,Ny,Jx,Jy,XDD,YDD);
    
            %Storing results
            FTLEF((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = FTLEf;
            
            % FTLE-backward
            [FTLEb] = ...
                FTLE_dgelem_f(XI,YI,-UI,-VI,XDGgt0,YDGgt0,dtftle,NTf,Nx,Ny,Jx,Jy,XDD,YDD);
            
            %Storing results
            FTLEB((j-1)*Ny+1:j*Ny+1,(i-1)*Nx+1:i*Nx+1) = FTLEb;
    
        end
    end
    MFTLEf(:,hh) = reshape(FTLEF,[],1);
    MFTLEb(:,hh) = reshape(FTLEB,[],1);
    hh = hh+1;
    end
    tt = tt+1;

%     
%     figure(2)
%     % Plotting U velocity
% %     ax1 = subplot(3,1,1);
%     surf(X,Y,Uplot);
%     hold on
%     plot3(XDGg,YDGg,10*ones(size(XDGg)),'*','Color','k')
% %     hold on
% %     streamslice(X,Y,Uplot,Vplot,'k');
% %     axis([xmin xmax ymin ymax -1 1.2]);
%     axis([x0dg-0.2 x0dg+lxdg+0.2 y0dg-0.2 y0dg+lydg+0.2 -10 10]);
%     view(0,90);
%     caxis([-0.4 0.8])
% %     colormap(ax1,jet)
%     colormap jet
%     shading interp
%     axis equal
%     colorbar
%     title(strcat('U, Flat Plate of ',string(lp),', Re=',string(round(Re)),', Time =',string(round(time,1)))); xlabel('X'); ylabel('Y');
%     drawnow
%     hold off
%         
% %     % Plotting V velocity
% %     ax2 = subplot(3,1,2);
% %     surf(X,Y,Vplot);
% % %     hold on
% % %     streamslice(X,Y,Uplot,Vplot,'k');
% %     axis([xmin xmax ymin ymax -1 1]);
% %     view(0,90);
% %     caxis([-0.4 0.4])
% %     colormap(ax2,jet)
% %     shading interp
% %     axis equal
% %     colorbar
% %     title('V'); xlabel('X'); ylabel('Y');
% %     drawnow
% % %     hold off
% %     
% %     % Plotting pressure
% %     ax3 = subplot(3,1,3);
% %     surf(X,Y,Pplot);
% % %     hold on
% % %     streamslice(X,Y,Uplot,Vplot,'k');
% %     axis([xmin xmax ymin ymax -1000 1000]);
% %     view(0,90);
% % %     caxis([-0.4 0.1])
% %     colormap(ax3,gray)
% %     shading interp
% %     axis equal
% %     colorbar
% %     title('P'); xlabel('X'); ylabel('Y');
% %     drawnow
% % %     hold off
    
end
% semilogy(t,normeu); axis([0 tf 1E-3 1E1]); grid on; hold on;
% subplot(2,1,1); title('U'); xlabel('X'); ylabel('Y'); colorbar; shading interp;
% subplot(2,1,2); title('V'); xlabel('X'); ylabel('Y'); colorbar; shading interp;
% toc
% save('field_results_41-101_281_79_Re300.mat','MU','MV','MP','dt','dx','dy','lp','n','m','Re','T','X','Y','xmax','xmin','ymax','ymin')
save('field_results_FTLE_40-60_251_91_Re200.mat','MU','MV','MP','dt','dx','dy','lp','n','m','Re','T','X','Y','xmax','xmin','ymax','ymin','XG','YG','Nx','Ny','Nxg','Nyg','MFTLEf','MFTLEb','sizeg')

