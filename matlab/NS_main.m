% clear; clc;
clear all
%%
%Solver Burgers 2D, with primitive variables scheme for the nonlinear
%terms.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical paramters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu = 0.01;          %Kinematic viscosity (ISU)
rho = 1000;         %Density of water, kg/m3 (ISU)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Space parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of nodes in each direction
n = 81;                      %Number of nodes in X direction
m = 81;                      %Number of nodes in Y direction

%Boundaries locations of the domain
xmin = 0.0;       %Left boundary
ymin = 0.0;       %Bottom boundary
xmax = 1.0;       %Top boundary
ymax = 1.0;       %Right boundary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Courant Number
CFL = 0.95;                                    %CFL (nonlinear parameter)

% Initial time
to = 0;                                     %Initial time

% Final time
tf = 2;                                  %Final time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CG configuration
%

%Tolerance for the residual of Conjugate Gradient Method
tolcg = 1E-8;

%Maximum number of iterations for CG solver (in difussion step)
mniterd = 100;

%
% Options for 'smallest' SVD solution
%

% Tolerance for finding minimum singular vector related
tolsing = 1E-12;

%Number of iterations for inverse power method
sing = 1;

% Maximum number of iterations for the solver, to find minimum singular
% vector related
mniterm = 1000;

%
% SOR Solver configuration for the laplacian
%

% Over-Relaxation coefficient (w=1 for Gauss-Seidel)
w = 13.4523570058092 * exp(-0.206450260650164 * (n*m).^-0.434163866503769) - 11.4497834449085;

%Tolerance for the residual of SOR
tolsor = 1E-4;

%Maximum number of iterations for SOR solver (in Poisson step)
mniters = 1000;

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

if tolsor > (min(dx,dy)*min(dx,dy))/2
    tolsor = (min(dx,dy)*min(dx,dy))/2;
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial condition for velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Uo = zeros(m*n,1);
Uo(upbound) = 1.0;
Uo(n*m) = 1.0;
Uo(n*m-n+1) = 1.0;
Vo = zeros(m*n,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Temporal parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = (min(dx,dy)*CFL/(max(abs(Uo))));    %Size of the time step
T = round((tf-to)/dt);                      %Number of time steps
t = to+dt:dt:(T*dt)+to;                        %Computational time

%Reynolds Number
Re = max(Uo)*max([(xmax-xmin) (ymax-ymin)])/nu;

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

% % % Computing matrix K for diffusion with Kronecker products
% % K11 = spdiags(ones(m-2,1),0,(m-2),(m-2));
% % K12 = spdiags([ax*ones(n-2,1) ap*ones(n-2,1) ax*ones(n-2,1)],[-1 0 1],(n-2),(n-2));
% % K21 = spdiags(ones(m-2,2),[-1 1],(m-2),(m-2));
% % K22 = spdiags(ay*ones(n-2,1),0,n-2,n-2);
% % K1 = kron(K11,K12);
% % K2 = kron(K21,K22);
% % K = K1 + K2;
% % %spy(K)

% Computing K with CSC storage
[K11v,K11r,K11c] = csc_diag(ones(m-2,1),0);
h1 = [ax*ones(n-2,1) ap*ones(n-2,1) ax*ones(n-2,1)];
d1 = [-1 0 1];
[K12v,K12r,K12c] = csc_diag(h1,d1);
h2 = [ones(m-3,1) ones(m-3,1)];
d2 = [-1 1];
[K21v,K21r,K21c] = csc_diag(h2,d2);
[K22v,K22r,K22c] = csc_diag(ay*ones(n-2,1),0);

[K1v,K1r,K1c] = csc_kron(K11v,K11r,K11c,K12v,K12r,K12c);
[K2v,K2r,K2c] = csc_kron(K21v,K21r,K21c,K22v,K22r,K22c);

[Kv,Kr,Kc] = csc_sum(K1v,K1r,K1c,K2v,K2r,K2c);
% csc_spy(Kv,Kr,Kc) 

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

%
% Computing matrix L for Poisson Equation with Kronecker products
% Standard Matlab

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

%
% Computing L through CSC storage
%
[L11v,L11r,L11c] = csc_diag(ones(m,1),0);
h1 = [alx*ones(n,1) alp*ones(n,1) alx*ones(n,1)];
h1(1,3) = alx*2;
h1(n-1,1) = alx*2;
d1 = [-1 0 1];
[L12v,L12r,L12c] = csc_diag(h1,d1);
h2 = ones(m-1,2);
h2(1,2) = 2;
h2(m-1,1) = 2;
d2 = [-1 1];
[L21v,L21r,L21c] = csc_diag(h2,d2);
[L22v,L22r,L22c] = csc_diag(aly*ones(n,1),0);

[L1v,L1r,L1c] = csc_kron(L11v,L11r,L11c,L12v,L12r,L12c);
[L2v,L2r,L2c] = csc_kron(L21v,L21r,L21c,L22v,L22r,L22c);

[Lv,Lr,Lc] = csc_sum(L1v,L1r,L1c,L2v,L2r,L2c);
%csc_spy(Lv,Lr,Lc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing minimum singular vector related of L' for regularization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % WITH STANDARD MATLAB FUNCTION
[UL1,SL1,VL1] = svds(L,1,'smallest');

% % % WITH CSC

% Calling preSOR
[Ltv,Ltr,Ltc] = csc_trans(Lv,Lr,Lc);
[LPv,LPr,LPc,LQv,LQr,LQc] = csc_preSOR(Lv,Lr,Lc,w);

%
% Computing smallest related singular vector with inverse power method
% C. POZRIKIDIS, 2001
% A Note on the Regularization of the Discrete Poisson–Neumann Problem
%

% Right Hand Side Singular Vector
VL = (1/sqrt(n*m))*ones(n*m,1);

% Left hand side singular vector (used for regularization)
% UL0 = 0.5*ones(n*m,1);
UL0 = -VL;
b = zeros(n*m,1);
nitersing = zeros(sing,1);
% normi = zeros(sing,1);

% LUs = zeros(length(Ltc)-1,1); pcs = 0;
% LUs = csc_diaga(Ltv,Ltr,Ltc); pcs = 1;
% LUs = abs(csc_diaga(Ltv,Ltr,Ltc)); pcs = 2;
% wsing = 1; LUs = csc_preconSSOR(Ltv,Ltr,Ltc,wsing); pcs = 3;
% % LUs = csc_preconILU0(Ltv,Ltr,Ltc); pcs = 4;
[Bv,Bc,Br] = css_trans(Ltv,Ltr,Ltc); pcs = 4;
[LUs] = csr_preconILU0(Bv,Bc,Br);
LUs = css_trans(LUs,Bc,Br);


r0as = b - csc_matvec(Ltv,Ltr,Ltc,UL0);
r0as = csc_matvec(Ltv,Ltr,Ltc,r0as);
kss = 8;
for i=1:sing
%     [UL,pit] = csc_CG(Ltv,Ltr,Ltc,b,UL0,mniterm,tolsing,LUs,pcs);
%     [UL,pit,ress] = csc_bicgstab(Ltv,Ltr,Ltc,b,UL0,r0as,mniterm,tolsing,LUs,pcs);
    [UL,pit,ress] = csc_gmres(Ltv,Ltr,Ltc,b,UL0,kss,mniterm,tolsing,LUs,pcs);

    UL = UL./norm(UL);
    nitersing(i) = pit;
    % check performance again matlab standard solver
    %  normi(i) = norm((abs(UL)-abs(UL1)))
end
pit
ress
norm(abs(UL)-abs(UL1))

%Regularization matrix
RM = (eye(n*m)-UL1*UL1');
rhs = ones(n*m,1);
rhs2 = RM*rhs;
rhtu = rhs'*UL;
rhs = rhs - rhtu*UL;
% norm(rhs-rhs2)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %
% % To store Angle not corrected (before Regularization)
% %
% angnc = zeros(length(t),1);
% 
% %
% % To store Angle corrected (after Regularization)
% %
% angc = zeros(length(t),1);
% 
% %
% % For storing SOR solver iterations taken
% %
% soriter = zeros(length(t),4);
% 
% j = 1;
% 
% %%
% %Temporal loop
% for time=t
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % First fractional step (Non-lineal advection)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %Derivatives
%     [dudx,dudy,dvdx,dvdy] = diver2D(Uo,Vo,dx,dy,n,m,bound,...
%         upboundin,doboundin,leboundin,riboundin);
%     
%     %Resultado (Explicito)
%     Up = Uo - dt*(Uo.*dudx + Vo.*dudy);
%     Vp = Vo - dt*(Uo.*dvdx + Vo.*dvdy);
%     
% %     %Boundary conditions
% %     Up(upbound) = 1.0;
% %     Up(n*m) = 1.0;
% %     Up(n*m-n+1) = 1.0;
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %  Second fractional step (Poisson equation)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %Derivacion de U y V resultado del primer paso fraccionado en sus respectivas direcciones principales
%     [dupdx,dvpdy] = grad2D(Up,Vp,dx,dy,n,m,upbound,dobound,lebound,ribound);
%     
%     %Vector de mano derecha sin regularizar
%     if(dy>=dx)
%         rhsp = ((dx^2)*rho/dt)*(dupdx + dvpdy);
%     else
%         rhsp = ((dy^2)*rho/dt)*(dupdx + dvpdy);
%     end
%         
%     %Deviation angle (just for test)
%     tetasr = acos(dot(rhsp,UL)*(1/(norm(rhsp)*norm(UL))));
%     angnc(j)=tetasr;
%     
%     %
%     %%% Regularizacion (Posrikidiz, 2001) %%%
%     %
%     
%     % Regularized right hand side
%     rhspr = RM*rhsp;
%     
%     % "Corrected angle" (just for test)
% %     tetar = acos(dot(rhspr,UL)*(1/(norm(rhspr)*norm(UL))));
%     tetar = dot(UL',rhspr);
%     angc(j) = tetar;
%         
%     %%% Neumann Boundary Condition == 0 for Poisson Equation %%%
%     %%% They are already included in matrix L structure
%     %%% Since Neumann boundary condition == 0, there's nothing to do with
%     %%% the right hand side (due to the use of central difference and ghost
%     %%% nodes for taking into account Neumann boundary condition)
%     
%     %
%     %%% Solving poisson equation %%%
%     %
%     
%     % % % With Matlab Standard function
%     % % P = L\rhspr;
%     
%     % With own csc_SOR solver
%     if j==1
%         P = zeros(n*m,1);
%     end
%     [P,nt] = csc_SOR(Lv,Lr,Lc,LPv,LPr,LPc,LQv,LQr,LQc,rhspr,P,mniters,tolsor);
%     soriter(j,1) = nt;
%     j=j+1;
%     
%     % Computing pressure gradient
%     [dpdx,dpdy] = grad2D(P,P,dx,dy,n,m,upbound,dobound,lebound,ribound);
%     Upp = Up - (dt/rho)*dpdx;
%     Vpp = Vp - (dt/rho)*dpdy;
% 
%     % % % Checking SOR performance in reference to Matlab Standard solver
%     % % [dpdx3,dpdy3] = grad2D(P3,P3,dx,dy,n,m,upbound,dobound,lebound,ribound);
%     % % soriter(j,3) = norm(abs(dpdx-dpdx3));
%     % % soriter(j,4) = norm(abs(dpdy-dpdy3));
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %  Third fractional step
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % Taking internal portion of Upp and Vpp (diffusion RHS) and settling
%     % dirichlet boundary conditions in the RHS
%     drhsx = zeros((n-2)*(m-2),1);
%     drhsy = zeros((n-2)*(m-2),1);
%     k=1;
%     for i=1:n*m
%         if ismember(i,bound)
%             continue
%         elseif (ismember(i,upboundin))
%             drhsx(k) = Upp(i) - ay*1;
%             drhsy(k) = Vpp(i);
%             k = k + 1;
%         else
%             drhsx(k) = Upp(i);
%             drhsy(k) = Vpp(i);
%             k = k + 1;
%         end
%     end
%     
% %     % Solving diffusion step with Standard Matlab solver
% %     Ud = K\drhsx;
% %     Vd = K\drhsy;
%     
%     % Solving diffusion step with Conjugate Gradient CSC method
%     [Ud,cgiu] = csc_CG(Kv,Kr,Kc,drhsx,drhsx,mniterd,tolcg);
%     [Vd,cgiv] = csc_CG(Kv,Kr,Kc,drhsy,drhsy,mniterd,tolcg);
%     
%     % fprintf('cgiu = %.0f \ncgiv = %.0f \n',cgiu,cgiv);
%         
%     % Adding boundary conditions
%     U = zeros(n*m,1);
%     V = zeros(n*m,1);
%     k = 1;
%     for i=1:n*m
%         if ismember(i,bound)
%             if (ismember(i,upbound) || i==(n*m-n+1) || i==n*m)
%                 U(i) = 1.0;
%             end
%         else
%             U(i) = Ud(k);
%             V(i) = Vd(k);
%             k = k + 1;
%         end
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Updating Uo and Vo
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Uo = U;
%     Vo = V;
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Postprocessing and plot
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     Uplot = reshape(Uo,n,m);
%     Uplot = Uplot';
%     Uplot = Uplot(m:-1:1,:);
%     Vplot = reshape(Vo,n,m);
%     Vplot = Vplot';
%     Vplot = Vplot(m:-1:1,:);
%         
%     % Plotting U velocity
%     subplot(1,2,1)
%     surf(X,Y,Uplot); 
%     axis([xmin xmax ymin ymax -1 1]);
%     view(0,90);
% %     caxis([-0.4 1])
%     shading interp
%     colorbar
%     drawnow
%         
%     % Plotting V velocity
%     subplot(1,2,2)
%     surf(X,Y,Vplot); 
%     axis([xmin xmax ymin ymax -1 1]);
%     view(0,90);
% %     caxis([-0.4 0.1])
%     shading interp
%     colorbar
%     drawnow
%     
% end
% % semilogy(t,normeu); axis([0 tf 1E-3 1E1]); grid on; hold on;
% subplot(1,2,1); title('U'); xlabel('X'); ylabel('Y'); colorbar; shading interp;
% subplot(1,2,2); title('V'); xlabel('X'); ylabel('Y'); colorbar; shading interp;



