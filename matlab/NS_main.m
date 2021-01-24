clear; clc;
%%
%Solver Burgers 2D, with primitive variables scheme for the nonlinear
%terms.

%%
%Entry constant values
nu = 0.0001;          %Kinematic viscosity 0.02 >= nu >= 0.0003
rho = 1000;         %Density of water, kg/m3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Space parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of nodes in each direction
n = 101;                     %Number of nodes in X direction
m = 101;                     %Number of nodes in Y direction

%Boundaries locations of the domain
xmin = 0.0;       %Left boundary
ymin = 0.0;       %Bottom boundary
xmax = 1.0;       %Top boundary
ymax = 1.0;       %Right boundary

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
Mpos2 = reshape(pos,n,m);
Mpos2 = Mpos2';
Mpos2 = Mpos2(m:-1:1,:);

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

%%
%Initial condition for velocity
Uo = zeros(m*n,1);
Uo(upbound) = 1.0;
Uo(n*m) = 1.0;
Uo(n*m-n+1) = 1.0;
Vo = zeros(m*n,1);

% % % Plotting initial condition
% % Uplot = reshape(Uo,n,m);
% % Uplot = Uplot';
% % Uplot = Uplot(m:-1:1,:);
% % Vplot = reshape(Vo,n,m);
% % Vplot = Vplot';
% % Vplot = Vplot(m:-1:1,:);
% % 
% % % Plotting U velocity
% % subplot(1,2,1)
% % surf(X,Y,Uplot); 
% % axis([xmin xmax ymin ymax -1 1]);
% % view(0,90);
% % drawnow
% % 
% % % Plotting V velocity
% % subplot(1,2,2)
% % surf(X,Y,Vplot); 
% % axis([xmin xmax ymin ymax -1 1]);
% % view(0,90);
% % drawnow

%%
%Temporal parameters
CFL = 1;                                    %CFL (nonlinear parameter)
to = 0;                                     %Initial time
tf = 2;                                  %Final time
dt = (min(dx,dy)*CFL/(max(abs(Uo))));    %Size of the time step
T = round((tf-to)/dt);                      %Number of time steps
t = to+dt:dt:(T*dt)+to;                        %Computational time

%Reynolds Number
Re = max(Uo)*max([(xmax-xmin) (ymax-ymin)])/nu;

%%
%Linear Parameters
Sx = ((dt*nu)/(dx^2));  %Linear for x
Sy = ((dt*nu)/(dy^2));  %Linear for y

%%
%Stiffnes linear difussion matrix
%constants
ap=(1+2*Sx+2*Sy);
ax=-Sx;
ay=-Sy;

% Diagonals for compute K with Kronecker products
K11 = spdiags(ones(m-2,1),0,(m-2),(m-2));
K12 = spdiags([ax*ones(n-2,1) ap*ones(n-2,1) ax*ones(n-2,1)],[-1 0 1],(n-2),(n-2));
K21 = spdiags(ones(m-2,2),[-1 1],(m-2),(m-2));
K22 = spdiags(ay*ones(n-2,1),0,n-2,n-2);
K1 = kron(K11,K12);
K2 = kron(K21,K22);
K = K1 + K2;
%spy(K)
 
%%
%Laplace operator matrix with five points second drivate and neuman
%boundary conditions

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
%spy(L)


% Autovecor asociado al autovalor = 0
% [vecprop0,valprop0] = svds(L,1,'smallest');
[UL,SL,VL] = svds(L,1,'smallest');
% V = (1/sqrt(n*m))*ones(n*m,1);


%%
%Temporal loop
cont = 1;
for time=t
    
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
        
    %Deviation angle (just for test)
    tetasr = acos(dot(rhsp,UL)*(1/(norm(rhsp)*norm(UL))));
    
    %
    %%% Regularizacion (Posrikidiz, 2001) %%%
    %
    
    % Regularized right hand side
    rhspr = (1-UL.^2).*rhsp;
    
    % Substracting the mean
    rhspr = rhspr - mean(rhspr);
    
    % "Corrected angle" (just for test)
    tetar = acos(dot(rhspr,UL)*(1/(norm(rhspr)*norm(UL))));
    
    %%% Neumann Boundary Condition == 0 for Poisson Equation %%%
    %%% They are already included in matrix L structure
    %%% Since Neumann boundary condition == 0, there's nothing to do with
    %%% the right hand side (due to the use of central difference and ghost
    %%% nodes for taking into account Neumann boundary condition)
    
    %%% Solving poisson equation %%%
    P = L\rhspr;
        
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
        if ismember(i,bound)
            continue
        elseif (ismember(i,upboundin))
            drhsx(k) = Upp(i) - ay*1;
            drhsy(k) = Vpp(i);
            k = k + 1;
        else
            drhsx(k) = Upp(i);
            drhsy(k) = Vpp(i);
            k = k + 1;
        end
    end
    
    % Solving diffusion step
    Ud = K\drhsx;
    Vd = K\drhsy;
    
    % Adding boundary conditions
    U = zeros(n*m,1);
    V = zeros(n*m,1);
    k = 1;
    for i=1:n*m
        if ismember(i,bound)
            if (ismember(i,upbound) || i==(n*m-n+1) || i==n*m)
                U(i) = 1.0;
            end
        else
            U(i) = Ud(k);
            V(i) = Vd(k);
            k = k + 1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Uo and Vo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Uo = U;
    Vo = V;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Postprocessing and plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Uplot = reshape(Uo,n,m);
    Uplot = Uplot';
    Uplot = Uplot(m:-1:1,:);
    Vplot = reshape(Vo,n,m);
    Vplot = Vplot';
    Vplot = Vplot(m:-1:1,:);
        
    % Plotting U velocity
    subplot(1,2,1)
    surf(X,Y,Uplot); 
    axis([xmin xmax ymin ymax -1 1]);
    view(0,90);
    drawnow
        
    % Plotting V velocity
    subplot(1,2,2)
    surf(X,Y,Vplot); 
    axis([xmin xmax ymin ymax -1 1]);
    view(0,90);
    drawnow
    
end
% semilogy(t,normeu); axis([0 tf 1E-3 1E1]); grid on; hold on;
subplot(1,2,1); title('U'); xlabel('X'); ylabel('Y'); colorbar; shading interp;
subplot(1,2,2); title('V'); xlabel('X'); ylabel('Y'); colorbar; shading interp;



