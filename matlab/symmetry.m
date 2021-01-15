% Proof of symmetry of matrices
clear; clc;
%%
%Solver Burgers 2D, with primitive variables scheme for the nonlinear
%terms.

%%
%Entry constant values
nu = 0.01;          %Kinematic viscosity 0.02 >= nu >= 0.0003
dens = 1000;         %Density of water, kg/m3

%%
n = 9;                     %Number of nodes in X direction
m = 9;                     %Number of nodes in Y direction
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
boundiup = Mpos(2,3:m-2);
boundidown = Mpos(n-1,3:m-2);
boundileft = Mpos(3:n-2,2)';
boundiright = Mpos(3:n-2,m-1)';
boundicor = [Mpos(2,2),Mpos(2,m-1),Mpos(n-1,2),Mpos(n-1,m-1)];
boundint = [Mpos(2,2:m-1), Mpos(n-1,2:m-1), Mpos(3:n-2,2)', Mpos(3:n-2,m-1)'];
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
K = zeros((n-2)*(m-2),(n-2)*(m-2));
iMpos = zeros(m,n);
ipos = 1:(n-2)*(m-2);
i2Mpos = reshape(ipos,(n-2),(m-2));
iMpos(2:m-1,2:n-1) = i2Mpos;

for i=1:m
    for j=1:n
        if ismember(Mpos(i,j),bound)
            continue;
        elseif ismember(Mpos(i,j),boundicor)
            K(iMpos(i,j),iMpos(i,j)) = -(1 + 2*Sx + 2*Sy);
            if iMpos(i,j-1) ~=0
               K(iMpos(i,j),iMpos(i,j-1)) = Sx;
            end
            if iMpos(i,j+1) ~=0
                K(iMpos(i,j),iMpos(i,j+1)) = Sx;
            end
            if iMpos(i-1,j) ~=0
                K(iMpos(i,j),iMpos(i-1,j)) = Sy;
            end
            if iMpos(i+1,j) ~=0
                K(iMpos(i,j),iMpos(i+1,j)) = Sy;
            end
        elseif ismember(Mpos(i,j),boundidown) || ismember(Mpos(i,j),boundiup)
            K(iMpos(i,j),iMpos(i,j)) = -(1 + 2*Sy + 30*Sx/12);
            if iMpos(i-1,j) ~=0
                K(iMpos(i,j),iMpos(i-1,j)) = Sy;
            end
            if iMpos(i+1,j) ~=0
                K(iMpos(i,j),iMpos(i+1,j)) = Sy;
            end
            if iMpos(i,j-1) ~=0
               K(iMpos(i,j),iMpos(i,j-1)) = (16*Sx/12);
            end
            if iMpos(i,j+1) ~=0
                K(iMpos(i,j),iMpos(i,j+1)) = (16*Sx/12);
            end
            if iMpos(i,j-2) ~=0
               K(iMpos(i,j),iMpos(i,j-2)) = -(Sx/12);
            end
            if iMpos(i,j+2) ~=0
                K(iMpos(i,j),iMpos(i,j+2)) = -(Sx/12);
            end
        elseif ismember(Mpos(i,j),boundileft) || ismember(Mpos(i,j),boundiright)
            K(iMpos(i,j),iMpos(i,j)) = -(1 + 30*Sy/12 + 2*Sx);
            if iMpos(i-1,j) ~=0
                K(iMpos(i,j),iMpos(i-1,j)) = (16*Sy/12);
            end
            if iMpos(i+1,j) ~=0
                K(iMpos(i,j),iMpos(i+1,j)) = (16*Sy/12);
            end
            if iMpos(i-2,j) ~=0
                K(iMpos(i,j),iMpos(i-2,j)) = -(Sy/12);
            end
            if iMpos(i+2,j) ~=0
                K(iMpos(i,j),iMpos(i+2,j)) = -(Sy/12);
            end
            if iMpos(i,j-1) ~=0
                K(iMpos(i,j),iMpos(i,j-1)) = Sx;
            end
            if iMpos(i,j+1) ~=0
                K(iMpos(i,j),iMpos(i,j+1)) = Sx;
            end
        else
            K(iMpos(i,j),iMpos(i,j)) = -(1 + 30*Sx/12 + 30*Sy/12);
            if iMpos(i,j-1) ~=0
               K(iMpos(i,j),iMpos(i,j-1)) = (16*Sy/12);
            end
            if iMpos(i,j+1) ~=0
                K(iMpos(i,j),iMpos(i,j+1)) = (16*Sy/12);
            end
            if iMpos(i-1,j) ~=0
                K(iMpos(i,j),iMpos(i-1,j)) = (16*Sx/12);
            end
            if iMpos(i+1,j) ~=0
                K(iMpos(i,j),iMpos(i+1,j)) = (16*Sx/12);
            end
            if iMpos(i,j-2) ~=0
               K(iMpos(i,j),iMpos(i,j-2)) = -(Sy/12);
            end
            if iMpos(i,j+2) ~=0
                K(iMpos(i,j),iMpos(i,j+2)) = -(Sy/12);
            end
            if iMpos(i-2,j) ~=0
                K(iMpos(i,j),iMpos(i-2,j)) = -(Sx/12);
            end
            if iMpos(i+2,j) ~=0
                K(iMpos(i,j),iMpos(i+2,j)) = -(Sx/12);
            end
        end
    end
end
            
H = K-K';
norm(H)
% No symmetry at all





