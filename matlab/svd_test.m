clear; clc;
%%
%Analysis of Singular vectors of matrix L

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Space parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of nodes in each direction
N = 100;                     %Number of nodes in X direction
M = 100;                     %Number of nodes in Y direction

% RES = zeros((N/10)*(M/10),5);    
% k=1;
% for i=10:10:N
%     for j=10:10:M
%         
%         n= i;
%         m = j;
%         
%         %Boundaries locations of the domain
%         xmin = 0.0;       %Left boundary
%         ymin = 0.0;       %Bottom boundary
%         xmax = 1.0;       %Top boundary
%         ymax = 1.0;       %Right boundary
% 
%         %Space differentials
%         dx = (xmax-xmin)/(n-1);           %X diferential
%         dy = (ymax-ymin)/(m-1);           %Y diferential
% 
%         %Laplace operator matrix with five points second drivate and neuman
%         %boundary conditions
% 
%         % Storing diagonals in one matrix and their positions
%         % Ensuring 
%         if(dy>=dx)
%             alx = 1;
%             aly = (dx^2)/(dy^2);
%             alp = -(2 + ((dx^2)/(dy^2))*2);
%         else
%             alx = (dy^2)/(dx^2);
%             aly = 1;
%             alp =-(2*((dy^2)/(dx^2)) + 2);
%         end
% 
%         L11 = spdiags(ones(m,1),0,m,m);
%         h1 = [alx*ones(n,1) alp*ones(n,1) alx*ones(n,1)];
%         h1(2,3) = alx*2;
%         h1(n-1,1) = alx*2;
%         L12 = spdiags(h1,[-1 0 1],n,n);
%         h2 = ones(m,2);
%         h2(2,2) = 2;
%         h2(m-1,1) = 2;
%         L21 = spdiags(h2,[-1 1],m,m);
%         L22 = spdiags(aly*ones(n,1),0,n,n);
%         L1 = kron(L11,L12);
%         L2 = kron(L21,L22);
%         L = L1 + L2;
%         %spy(L)
% 
% 
%         % Autovecor asociado al autovalor = 0
%         % [vecprop0,valprop0] = svds(L,1,'smallest');
%         [U1,S1,V1] = svds(L,1,'smallest');
%         RES(k,1) = n;
%         RES(k,2) = m;
%         RES(k,3) = max(abs(U1));
%         RES(k,4) = S1;
%         RES(k,5) = max(abs(V1));
%         k = k + 1;
%     end
% end

load('svms.mat')

MS = reshape(RES(:,1),(N/10),(M/10));
NS = reshape(RES(:,2),(N/10),(M/10));
US = reshape(RES(:,3),(N/10),(M/10));
VS = reshape(RES(:,5),(N/10),(M/10));

surf(log(NS),log(MS),log(US))
xlabel('N')
ylabel('M')
zlabel('U')

% hold on
% 
% surf(log(NS),log(MS),log(VS))
% xlabel('N')
% ylabel('M')
% zlabel('V')
% 
% surf(log(NS),log(MS),log(US.*(MS.*NS)))
% xlabel('N')
% ylabel('M')
% zlabel('U')
% % 
% surf(log(NS),log(MS),log(VS.*(MS.*NS)))
% xlabel('N')
% ylabel('M')
% zlabel('V')

UMN = (US.*(MS.*NS));
LUS = log(US);
LNS = log(NS);
LMS = log(MS);
m = (LUS(10,1)-LUS(1,1))/(LNS(10,1)-LNS(1,1));
b = (LUS(10,1))-m*(LNS(10,1));

US2 = (NS.^m).*(MS.^m);

plot(LNS(:,1),LUS(:,1))
hold on
plot(LMS(1,:),LUS(1,:))
xlabel('log(dim)')
ylabel('log(U)')

VMN = (VS.*(MS.*NS)).^2;
VS2 = 1./(sqrt(MS.*NS));
