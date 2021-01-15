clear; clc;
%%
%Solver Burgers 2D, with primitive variables scheme for the nonlinear
%terms.

%%
%Entry constant values
% cont2 = 1;
% for nu = [0.02,0.008,0.002,0.0003]
nu = 0.02;          %Kinematic viscosity 0.02 >= nu >= 0.0003 

%%
% cont2 = 1;
% for num = [21, 51, 81, 101, 151]
%Spatial parameters
n = 101;                     %Number of nodes in X direction
m = 101;                     %Number of nodes in Y direction
% n = num;
% m = num;
a = -1;                     %Left and down boundary conditions
b = 1;                      %Rigth and up boundary conditions
dx = (b-a)/(n-1);           %X diferential
dy = (b-a)/(m-1);           %Y diferential
x = -1:dx:1;                %x vector
y = -1:dy:1;                %y vector
invy = 1:-dy:-1;            %y vector turned
[X,Y] = meshgrid(x,invy);   %Coord. matrices for x and y vectors
pos = 1:m*n;
Mpos = reshape(pos,n,m);

%Boundary nodes
bound = [Mpos(:,1)',Mpos(:,m)',Mpos(1,2:m-1),Mpos(n,2:m-1)];
bound = sort(bound);
boundint = [Mpos(2,2:m-1), Mpos(n-1,2:m-1), Mpos(3:n-2,2)',Mpos(3:n-2,m-1)'];
boundint = sort(boundint);
boundup = Mpos(1,:);
bounddown = Mpos(n,:);
boundleft = Mpos(:,1)';
boundright = Mpos(:,m)';

%%
%Initial condition for velocity
Uo = -sin(pi*X);
Vo = -(-sin(pi*Y));

%%
%Temporal parameters
CFL = 0.25;                                    %CFL (nonlinear parameter)
to = 0;                                     %Initial time
tf = 2/pi;                                  %Final time
dt = (min(dx,dy)*CFL/(max(max(abs(Uo)))));    %Size of the time step
T = round((tf-to)/dt);                      %Number of time steps
t = to+dt:dt:(T*dt)+to;                        %Computational time

%%
%Parameters
Sx = ((dt*nu)/(dx^2));  %Linear for x
Sy = ((dt*nu)/(dy^2));  %Linear for y

%%
%Stiffnes linear difussion matrix
diag0k = -(1 + 30*Sx/12 + 30*Sy/12)*ones(1,n*m);
diag1k = (16*Sy/12)*ones(1,n*m-1);
diagmk = (16*Sx/12)*ones(1,n*m-m);
diag2k = -(Sy/12)*ones(1,n*m-2);
diag2mk = -(Sx/12)*ones(1,n*m-2*m);
K = sparse(diag(diag0k,0)) + sparse(diag(diag1k,1)) + sparse(diag(diag1k,-1)) + sparse(diag(diagmk,m)) + sparse(diag(diagmk,-m) + sparse(diag(diag2k,2)) + sparse(diag(diag2k,-2)) + sparse(diag(diag2mk,2*m)) + sparse(diag(diag2mk,-2*m)));
% K = sparse(K);
K(bound,:) = 0;
K(boundint,:) = 0;

for i=boundint
    K(i,i) = -(1 + 2*Sx + 2*Sy);
    K(i,i-m) = Sx;
    K(i,i+m) = Sx;
    K(i,i-1) = Sy;
    K(i,i+1) = Sy;
end
Ku = K;
Kv = K;
for i=boundup
    Kv(i,i) = 1;
    Ku(i,i) = 1;
%     Ku(i,i) = 1/dy;
%     Ku(i,i+1) = -1/dy;
end
for i=bounddown
    Kv(i,i) = 1;
    Ku(i,i) = 1;
%     Ku(i,i) = -1/dy;
%     Ku(i,i-1) = 1/dy;
end
for i=boundleft
    Ku(i,i) = 1;
    Kv(i,i) = 1;
%     Kv(i,i) = -1/dx;
%     Kv(i,i+n) = 1/dx;
end
for i=boundright
    Ku(i,i) = 1;
    Kv(i,i) = 1;
%     Kv(i,i) = 1/dy;
%     Kv(i,i-n) = -1/dy;
end
%%
%Temporal loop
cont = 1;
for k=t
    %Anlitica
    [Ua,Va] = analitica2D(X,Y,k,nu);
    Va = -Va;
    %Primer paso fraccionado, adveccion explicita
    [DMUp1,DMUp2,DMVp1,DMVp2] = derivada2D(Uo,Vo,dx,dy);
    MUp = Uo - dt*Uo.*DMUp1 - dt*Vo.*DMUp2;
    MVp = Vo - dt*Uo.*DMVp1 - dt*Vo.*DMVp2;
%     [MUp,MVp] = deriv2D_up2(Uo,Vo,X,Y,dt);
    %Condiciones de frontera para paso fraccionado 2
    MUp(bound) = -Ua(bound);
    MVp(bound) = -Va(bound);
    %Vectorizacion de las matrices obtenidas en paso fraccionado 1
    Up = MUp(:);
    Vp = MVp(:);
    %Segundo paso fraccionado
    vecU = Ku\-(Up);
    vecV = Kv\-(Vp);
    %Reshape soluciones en matrices
    U = reshape(vecU,m,n);
    V = reshape(vecV,m,n);
    %Postproceso
    difeu = U-Ua;
    difev = V-Va;
    normeu(cont) = norm(difeu,Inf);
    normev(cont) = norm(difev,Inf);
    cont = cont+1;
    Uo = U;
    Vo = V;
    %subplot(2,2,1)
%     surf(X,Y,(abs(difeu))); axis([-1 1 -1 1 1E-6 1])
%     drawnow
    %grid on
%     subplot(2,2,2)
%     surf(X,Y,V); axis([-1 1 -1 1 -2 2]); %hold on;
%     drawnow
%     surf(X,Y,Va); 
%     drawnow; hold off
%     hold on;
%     surf(X,Y,Ua); axis([-1 1 -1 1 -2 2]); hold off;
%     drawnow
%     %grid on
%     subplot(2,2,3)
%     surf(X,Y,(abs(difev))); axis([-1 1 -1 1 1E-6 1])
%     drawnow
%     %grid on
%     subplot(2,2,4)
%     surf(X,Y,V); axis([-1 1 -1 1 -2 2])
%     drawnow
%     %grid on
%     %pause(dt)
end
semilogy(t,normeu); axis([0 tf 1E-3 1E1]); grid on; hold on;
% numberit(1,cont2)=n*m;
% numberit(2,cont2)=dt;
% numberit(3,cont2)=normeu(T);
% numberit(4,cont2)=normev(T);
% numberit(5,cont2)=nu;
% cont2=cont2+1;
% end