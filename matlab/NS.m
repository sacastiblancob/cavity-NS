clear; clc;
%%
%Solver Burgers 2D, with primitive variables scheme for the nonlinear
%terms.

%%
%Entry constant values
nu = 0.01;          %Kinematic viscosity 0.02 >= nu >= 0.0003
dens = 1000;         %Density of water, kg/m3

%%
%Spatial parameters
n = 201;                     %Number of nodes in X direction
m = 201;                     %Number of nodes in Y direction
% n = num;
% m = num;
a = 0;                     %Left and down boundary conditions
b = 1;                      %Rigth and up boundary conditions
dx = (b-a)/(n-1);           %X diferential
dy = (b-a)/(m-1);           %Y diferential
x = a:dx:b;                %x vector
y = a:dy:b;                %y vector
invy = b:-dy:a;            %y vector turned
[X,Y] = meshgrid(x,invy);   %Coord. matrices for x and y vectors
pos = 1:m*n;
Mpos = reshape(pos,n,m);

%Boundary nodes
bound = [Mpos(:,1)',Mpos(:,m)',Mpos(1,2:m-1),Mpos(n,2:m-1)];
bound = sort(bound);
boundup = Mpos(1,1:m);
bounddown = Mpos(n,1:m);
boundleft = Mpos(2:n-1,1)';
boundright = Mpos(2:n-1,m)';
boundint = [Mpos(2,2:m-1), Mpos(n-1,2:m-1), Mpos(3:n-2,2)',Mpos(3:n-2,m-1)'];
boundint = sort(boundint);

%%
%Initial condition for velocity
Uo = sparse(n,m);
Uo(1,:) = 1;
Vo = sparse(n,m);
% Uo = -sin(pi*X);
% Vo = -(-sin(pi*Y));

%%
%Temporal parameters
CFL = 1;                                    %CFL (nonlinear parameter)
to = 0;                                     %Initial time
tf = 2;                                  %Final time
dt = (min(dx,dy)*CFL/(max(max(abs(Uo)))));    %Size of the time step
T = round((tf-to)/dt);                      %Number of time steps
t = to+dt:dt:(T*dt)+to;                        %Computational time

%%
%Linear Parameters
Sx = ((dt*nu)/(dx^2));  %Linear for x
Sy = ((dt*nu)/(dy^2));  %Linear for y

%%
%Stiffnes linear difussion matrix

% Storing diagonals in one matrix and their positions
e = ones(n*m,1);
diags = [-(Sx/12)*e, (16*Sx/12)*e, -(Sy/12)*e, (16*Sy/12)*e, -(1 + 30*Sx/12 + 30*Sy/12)*e, ...
    (16*Sy/12)*e, -(Sy/12)*e, (16*Sx/12)*e, -(Sx/12)*e];
diagpos = [-2*m,-m,-2,-1,0,1,2,m,2*m];

% Making K
K = spdiags(diags,diagpos,m*n,m*n);
 
% Cleaning boundaries
K(bound,:) = 0;
K(boundint,:) = 0;

% Changing internal boundaries
for i=boundint
    K(i,i) = -(1 + 2*Sx + 2*Sy);
    K(i,i-m) = Sx;
    K(i,i+m) = Sx;
    K(i,i-1) = Sy;
    K(i,i+1) = Sy;
end
% Ku = K;
% Kv = K;

% Changing top boundaries
for i=boundup
    K(i,i) = 1;
%     Ku(i,i) = 1;
%     Ku(i,i) = 1/dy;
%     Ku(i,i+1) = -1/dy;
end

% Changing bottom boundaries
for i=bounddown
    K(i,i) = 1;
%     Ku(i,i) = 1;
%     Ku(i,i) = -1/dy;
%     Ku(i,i-1) = 1/dy;
end

% Changing left boundaries
for i=boundleft
    K(i,i) = 1;
%     Kv(i,i) = 1;
%     Kv(i,i) = -1/dx;
%     Kv(i,i+n) = 1/dx;
end

% Changing right boundaries
for i=boundright
    K(i,i) = 1;
%     Kv(i,i) = 1;
%     Kv(i,i) = 1/dy;
%     Kv(i,i-n) = -1/dy;
end

%%
%Laplace operator matrix with five points second drivate and neuman
%boundary conditions

% Storing diagonals in one matrix and their positions
diags = [-(1/(12*dx^2))*e, (16/(12*dx^2))*e, -(1/(12*dy^2))*e, (16/(12*dy^2))*e, ...
    -30*(1/(12*dx^2) + 1/(12*dy^2))*e, (16/(12*dy^2))*e,...
    -(1/(12*dy^2))*e, (16/(12*dx^2))*e, -(1/(12*dx^2))*e];

% Computing laplacian matrix
L = spdiags(diags,diagpos,m*n,m*n);

% Cleaning boundaries
L(bound,:) = 0;
L(boundint,:) = 0;
for i=boundint
    L(i,i) = (-2/(dx^2)) + (-2/(dy^2));
    L(i,i-m) = 1/(dx^2);
    L(i,i+m) = 1/(dx^2);
    L(i,i-1) = 1/(dy^2);
    L(i,i+1) = 1/(dy^2);
end
for i=boundup
%     L(i,i) = 1/(dy);
%     L(i,i+1) = -1/(dy);
    L(i,i) = 3/(2*dy);
    L(i,i+1) = -4/(2*dy);
    L(i,i+2) = 1/(2*dy);
end
for i=bounddown
%     L(i,i) = -1/(dy);
%     L(i,i-1) = 1/(dy);
    L(i,i) = -3/(2*dy);
    L(i,i-1) = 4/(2*dy);
    L(i,i-2) = -1/(2*dy);
end
for i=boundleft
    L(i,i+2*m) = -1/(2*dx);
    L(i,i+m) = 4/(2*dx);
    L(i,i) = -3/(2*dx);
%     L(i,i) = -1/dx;
%     L(i,i+m) = 1/dx;
%     L(i,i+1) = 1/(2*dy);
%     L(i,i-1) = -1/(2*dy);
end
for i=boundright
    L(i,i-2*m) = 1/(2*dx);
    L(i,i-m) = -4/(2*dx);
    L(i,i) = 3/(2*dx);
%     L(i,i) = 1/dx;
%     L(i,i-m) = -1/dx;
%     L(i,i+1) = 1/(2*dy);
%     L(i,i-1) = -1/(2*dy);
end
% Autovecor asociado al autovalor = 0
[vecprop0,valprop0] = svds(L,1,'smallest');

% LF = zeros(n*m,n*m);
% for i=1:n*m
%     LF(i,:) = L(i,:);
% end

LF = L*L';

H = elgauss_sing(LF,zeros(n*m,1));

%%
%Temporal loop
cont = 1;
for k=t
    %
    %%%%%%  Primer paso fraccionado (Advección no lineal)   %%%%%%
    %
    %Derivadas
    [DMUpx,DMVpy,DMUpy,DMVpx] = derivada2D(Uo,Vo,dx,dy);
    
    %Resultado (Explicito)
    MUp = Uo - dt*Uo.*DMUpx - dt*Vo.*DMUpy;
    MVp = Vo - dt*Uo.*DMVpx - dt*Vo.*DMVpy;
    
%     MUp(bound) = -Ua(bound);
%     MVp(bound) = -Va(bound);
%     pause

    %
    %%%%%%  Segundo paso fraccionado (Poisson equation) %%%%%%
    %
    
    %Derivacion de U y V resultado del primer paso fraccionado en sus respectivas direcciones principales
    [DMUpx,DMVpy,DMUpy,DMVpx] = derivada2D(MUp,MVp,dx,dy);
    
    %Vector de mano derecha sin regularizar
    vecb = (dens/dt)*(DMUpx(:) + DMVpy(:));            %Vector de mano derecha sin regularizar de la ecuacion de Poisson
    
    %Angulos de desviación
    tetasr = acos(dot(vecb,vecprop0)*(1/(norm(vecb)*norm(vecprop0))));
    
    %%% Regularizacion %%%
    
    %Vector de mano derecha regularizado
    vecbr = (1-vecprop0.^2).*vecb;
    
    %Ángulos corregidos
    tetar = acos(dot(vecbr,vecprop0)*(1/(norm(vecbr)*norm(vecprop0))));
    
    %%% Condiciones de frontera neuman = 0 para la presion %%%
    vecbr(bound) = 0;
    
    %%% Solucionando la ecuacion de Poisson %%%
    % Presión vector
    P = L\vecbr;
        
    % Presión matriz
    MP = reshape(P,n,m);
    
    % Amarrando el resultado de la presión a algún valor determinado
    amarre = MP(1,1) - 1;           %Amarrando P a 1 en la esquina superior izquierda
    MPamarrada = MP - amarre;       %Amarrando toda la matriz
    
    % Cálculo del gradiente de la presión
    [DMPx,DMPy] = derivada2D(MP,MP,dx,dy);
    MUpp = MUp - (dt/dens)*DMPx;
    MVpp = MVp - (dt/dens)*DMPy;
    
%     pause
    %
    %%%%%%  Tercer paso fraccionado   %%%%%%
    %
    %Vectorizacion resultados de paso fraccionado 2
    Upp = MUpp(:);
    Vpp = MVpp(:);
    
    % Condiciones de frontera
    Upp(bound) = 0;
    Upp(boundup) = -1;
    Vpp(bound) = 0;
    
    % Solucion sistema lineal de ecuaciones de difusion
    vecU = K\-(Upp);
    vecV = K\-(Vpp);
    U = reshape(vecU,n,m);
    V = reshape(vecV,n,m);    

%     pause
    %%%%%%  Reasignacion   %%%%%%
    Uo = U;
    Vo = V;

    %%%%%%  Postproceso   %%%%%%
    V = -V;         %Esta sustitucion es necesaria para ver los resultados 
                    %pues en el calculo de la velocidad en Y, internamente es positiva hacia abajo
                    
%     [UV] = stream2(X,Y,U,V,0.55*ones(1,21),0:0.05:1,0.01);
    
    % Graficando solución velocidad U
    subplot(1,2,1)
    surf(X,Y,U); axis([a b a b -2 2]); view(0,90); %hold on;
    drawnow
    
    % Graficando solución velocidad V
    subplot(1,2,2)
    surf(X,Y,V); axis([a b a b -2 2]); view(0,90);%hold on;
    drawnow
    
    % Graficando solución presión amarrada
%     subplot(2,2,3)
%     streamline(UV); axis([a b a b]);
%     drawnow
%     surf(X,Y,MPamarrada); axis([a b a b 0 8000]); %hold on;
%     drawnow

end
% semilogy(t,normeu); axis([0 tf 1E-3 1E1]); grid on; hold on;
