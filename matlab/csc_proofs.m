% % B = [2 0 0 -5 4; 2 2 0 0 5; 0 0 6 0 1; 2 0 0 8 2; 0 2 0 1 5];
% A = [2 0 1 5 4; 2 2 0 0 5;0 0 6 1 1; 0 0 0 2 0;5 0 0 0 5];
% % %B = [1 6 9 12 15;2 7 0 0 16;3 0 10 0 17;4 0 0 13 18;5 8 11 14 19];
% % %  A = diag([1:5]);
% % %  A(5,1) = 6;
% % %  A(1,3) = 3;
% % %  A(4,2) = 4;
% % %  A(4,3) = 2;
% % %A = [2 0 1 5 4; 2 2 0 0 4;2 0 6 5 4; 0 0 0 2 0;5 0 0 0 5];
% % % n = length(A);
% b = [1;1;1;1;1];
% % % Aorg = A;
% % % y = b;
% % % x = zeros(n,1);
% % 
% [Av,Ar,Ac] = full2csc(A);
% % [Bv,Br,Bc] = full2csc(B);
% %
% 
% cr = A*b;
% 
% nz = length(Av);
% n = max(Ar);
% m = length(Ac)-1;
% % if m~=length(b)
% %     disp('ERROR!!! Dimensions does not agree')
% %     return
% % end
% c = zeros(n,1);
% 
% for j=1:m
%     for i=Ac(j):Ac(j+1)-1
%        r = Ar(i);
%        c(r) = c(r) + Av(i)*b(j);
%     end
% end
%        
        
        
% AB = A+B;
% [ABvo,ABro,ABco] = full2csc(AB);
% 
% [ABv,ABr,ABc] = csc_sum(Av,Ar,Ac,Bv,Br,Bc);
% 
% % load('Kmat5000')
% % [Kv,Kr,Kc] = full2csc(K);
% % tic
% % [ABv,ABr,ABc] = csc_sum(Kv,Kr,Kc,Kv,Kr,Kc);
% % toc

% d0 = [1;2;3;4;5];
% d1 = [6;7;8;9;10];
% dm1 = [6;7;8;9;10];
% h2 = [dm1 d0 d1];
% d = [-2 2 3];
% dp = d(abs(d) == min(abs(d)));
% dp = dp(1);
% m = length(h2(:,d==dp))+abs(dp);
% nz=0;
% for i=1:length(d)
%     %nz = nz+length(h2)-abs(d(i))+(m-length(h2));
%     nz = nz - abs(d(i)) + m;
% end
% Av = zeros(nz,1);
% Ar = zeros(nz,1); 
% Ac = zeros(m+1,1);
% Ac(1) = 1;
% k = 0;
% c = -m;
% cv = d*0;
% p = 1;
% for i=1:m
%     Ac(i+1) = Ac(i) + sum((d<=k).*(d>c));
%     Ar(Ac(i):Ac(i+1)-1) = sort(abs(d(logical((d<=k).*(d>c)))-1-k));
%     cv = cv+(d<=k).*(d>c);
%     hcols = h2(:,logical((d<=k).*(d>c)));
%     hrows = flip(cv(logical((d<=k).*(d>c))));
%     row = 1;
%     col = length(hrows);
%     for j=Ac(i):Ac(i+1)-1
%         Av(p) = hcols(hrows(row),col);
%         row = row + 1;
%         col = col - 1;
%         p = p + 1;
%     end
%     k=k+1;
%     c=c+1;
% end

% d0 = [1;2;3;4;5];
% d1 = [6;7;8;9;10];
% dm1 = [6;7;8;9;10];
% h2 = [dm1 d0 d1];
% d = [-2 2 3];
% 
% [Av,Ar,Ac] = csc_diag(h2,d);
% 
% d3 = ones(10,1);
% d0 = 2*ones(10,1);
% h3 = [d3 d3 d0 d3 d3];
% dia = [-2 -1 0 2 3];
% 
% [Bv,Br,Bc] = csc_diag(h3,dia);
% 
% csc_spy(Bv,Br,Bc)

% h3m = [h3;zeros(1,5)];
% B = spdiags(h3m,dia,11,11);
% [Bvo,Bro,Bco] = full2csc(B);
% norm(abs(Bv-Bvo))
% norm(abs(Br-Bro))
% norm(abs(Bc-Bco))

% A = [1 3 0 0; 2 4 0 0; 0 0 5 0; 0 0 0 6];
% B = [1 0 0; 1 2 1; 0 0 3];
% C = kron(A,B);
% [Av,Ar,Ac] = full2csc(A);
% [Bv,Br,Bc] = full2csc(B);
% [Cvo,Cro,Cco] = full2csc(C);
% [Cv,Cr,Cc] = csc_kron(Av,Ar,Ac,Bv,Br,Bc);

% ma = length(Ac)-1;
% mb = length(Bc)-1;
% nza = Ac(ma+1)-1;
% nzb = Bc(mb+1)-1;
% nb = max(Br);
% mc = ma*mb;
% nzc = nza*nzb;
% Cc = zeros(mc+1,1);
% Cv = zeros(nzc,1);
% Cr = zeros(nzc,1);
% Cc(1) = 1;
% aj=1;
% bj=1;
% p=1;
% for j=1:mc
%     af = (Ac(aj+1) - Ac(aj));
%     bf = (Bc(bj+1) - Bc(bj));
%     Cc(j+1) = Cc(j) + bf*af;
%     r = nb;
%     for ia = Ac(aj):Ac(aj+1)-1
%         %    for i=Cc(j):Cc(j+1)-1
%         for ib = Bc(bj):Bc(bj+1)-1
%             Cr(p) = r*(Ar(ia)-1)+Br(ib);
%             Cv(p) = Av(ia)*Bv(ib);
%             p = p+1;
%         end
%     end
%     bj=bj+1;
%     if bj==mb+1
%         aj=aj+1;
%         bj=1;
%     end
% end
% 
% %Building difussion-symmetric matrix
% 
% Dx=0.25;
% Dy=0.25;
% u=0;
% v=0;
% 
% %Discretizacion espacial
% a=4;
% b=4;
% dx=(8/72);
% dy=(8/72);
% 
% n=round((a-(-b))/dx+1);
% m=round((a-(-b))/dy+1);
% x=-a:dx:b;
% y=-a:dy:b;
% X=kron(ones(1,m),x);
% Y1=kron(y,ones(1,n));
% Y=Y1(n*m:-1:1);
% MX=reshape(X,m,n)';
% MY=reshape(Y,m,n)';
% 
% %Discretizacion temporal
% dt=.025;
% tinicial=1+dt;
% tfinal=10;
% t=tinicial+dt:dt:tfinal;
% 
% %Constantes de la ecuacion diferencial computacional
% Sx=(Dx*dt)/(dx^2);
% Sy=(Dy*dt)/(dy^2);
% CFLx=(u*dt)/dx;
% CFLy=(v*dt)/dy;
% CFL = max(CFLx,CFLy);

% %Valores de los nodos que estan en la frontera
% front=[1:n,n+1:n:n*m-2*n+1,2*n:n:n*m-n,n*m-n+1:n*m];
% front=sort(front);
% fronts=(n+3:2*n-2);
% frontd=(3*n-1:n:n*m-2*n-1);
% frontin=(n*m-2*n+3:n*m-n-2);
% frontiz=(2*n+2:n:n*m-3*n+2);

% %Ensamblaje matriz sistema de ecuaciones para el paso 2 del planteamiento
% %numerico
% ap=(1+2*Sx+2*Sy);
% ax=-Sx;
% ay=-Sy;
% 
% K11=diag(ones(n-2,1));
% K12=ax*diag(ones(n-3,1),1)+ax*diag(ones(n-3,1),-1)+ap*diag(ones(n-2,1));
% K21=diag(ones(n-3,1),1)+diag(ones(n-3,1),-1);
% K22=ay*diag(ones(n-2,1));
% K1=kron(K11,K12);
% K2=kron(K21,K22);
% K = K1 + K2;
% 
% [K11v,K11r,K11c] = csc_diag(ones(n-2,1),0);
% h1 = [ax*ones(n-2,1) ap*ones(n-2,1) ax*ones(n-2,1)];
% d1 = [-1 0 1];
% [K12v,K12r,K12c] = csc_diag(h1,d1);
% h2 = [ones(n-3,1) ones(n-3,1)];
% d2 = [-1 1];
% [K21v,K21r,K21c] = csc_diag(h2,d2);
% [K22v,K22r,K22c] = csc_diag(ay*ones(n-2,1),0);
% 
% [K1v,K1r,K1c] = csc_kron(K11v,K11r,K11c,K12v,K12r,K12c);
% [K2v,K2r,K2c] = csc_kron(K21v,K21r,K21c,K22v,K22r,K22c);
% 
% [Kv,Kr,Kc] = csc_sum(K1v,K1r,K1c,K2v,K2r,K2c);
% clearvars K11v K11r K11c K12v K12r K12c K21v K21r K21c K22v K22r K22c K1v K1r K1c K2v K2r K2c
% 
% caux = ones(length(Kc)-1,1);
% % A = [6 0 1 2 1; 0 4 0 0 -2;1 0 6 -1 0; 2 0 -1 5 0;1 -2 0 0 5];
% % b = [1;1;1;1;1];
% % %A = [2 1 0; 1 3 1; 0 1 1];
% % %b = [1;2;3];
% % [Av,Ar,Ac] = full2csc(A);
% x = caux;
% xo = x;
% c = K\caux;

% %Jacobi
% [Pv,Qv,Qc,Qr] = csc_prejacobi(Kv,Kr,Kc);
% [x,t] = csc_jacobi(Kv,Kr,Kc,Pv,Qv,Qc,Qr,caux,x,100,1E-4);

% %SOR
% [Pv,Pr,Pc,Qv,Qr,Qc] = csc_preSOR(Kv,Kc,Kr,1);
% [x,t] = csc_SOR(Kv,Kr,Kc,Pv,Pr,Pc,Qv,Qr,Qc,caux,x,100,1E-4);

% % CG
% [x,t] = csc_CG(Kv,Kr,Kc,caux,x,100,1E-4);

% % DPCG
% [x,t] = csc_DPCG(Kv,Kr,Kc,b,caux,100,1E-4);

% %      A x = b
% % P x(t+1) = (P - A) x(t) + b
% % P x(t+1) = Q x(t) + b   ;   Q = P - A
% 
% siz = size(A);
% n = siz(1);
% m = siz(2);
% x = b;
% L = zeros(n,m);
% U = zeros(n,m);
% D = zeros(n,m);
% for j=1:m
%     for i=1:n
%         if i<j
%             U(i,j) = A(i,j);
%         elseif i>j
%             L(i,j) = A(i,j);
%         else
%             D(i,j) = A(i,j);
%         end
%     end
% end


% %LUD CSC
% m = length(Ac)-1;
% n = max(Ar);
% nzu = 0;
% nzl = 0;
% nzd = 0;
% for j=1:m
%     for i=Ac(j):Ac(j+1)-1
%         if Ar(i)<j
%             nzu = nzu + 1;
%         elseif Ar(i)>j
%             nzl = nzl + 1;
%         else
%             nzd = nzd + 1;
%         end
%     end
% end

% Uv = zeros(nzu,1);
% Ur = zeros(nzu,1);
% Uc = zeros(m+1,1);
% Uc(1) = 1;
% Lv = zeros(nzl,1);
% Lr = zeros(nzl,1);
% Lc = zeros(m+1,1);
% Lc(1) = 1;
% Dv = zeros(nzd,1);
% Dr = zeros(nzd,1);
% Dc = zeros(m+1,1);
% Dc(1) = 1;
% up=1;
% lp=1;
% dp=1;
% for j=1:m
%     for i=Ac(j):Ac(j+1)-1
%         if Ar(i)<j
%             Uv(up) = Av(i);
%             Ur(up) = Ar(i);
%             up = up + 1;
%         elseif Ar(i)>j
%             Lv(lp) = Av(i);
%             Lr(lp) = Ar(i);
%             lp = lp + 1; 
%         else
%             Dv(dp) = Av(i);
%             Dr(dp) = Ar(i);
%             dp = dp + 1;
%         end
%     end
%     Uc(j+1) = up;
%     Lc(j+1) = lp;
%     Dc(j+1) = dp;
% end

% %Jacobi CSC
% 
% tol = 1E-4;
% 
% %Matrices
% m = length(Ac)-1;
% nzd = length(Av)-m;
% Qv = zeros(nzd,1);
% Qr = zeros(nzd,1);
% Qc = zeros(m+1,1);
% Qc(1) = 1;
% Pv = zeros(m,1);
% qp=1;
% pp=1;
% for j=1:m
%     for i=Ac(j):Ac(j+1)-1
%         if Ar(i) ~= j
%             Qv(qp) = -Av(i);
%             Qr(qp) = Ar(i);
%             qp = qp + 1;
%         else
%             Pv(pp) = Av(i);
%             pp = pp + 1;
%         end
%     end
%     Qc(j+1) = qp;
% end
% 
% %loop
% for t=1:100
%     d = csc_matvec(Qv,Qr,Qc,x) + b;
%     x = d./Pv;
%     if norm(b - csc_matvec(Av,Ar,Ac,x))<tol
%         return
%     end
% end

% %Jacobi
% tol = 1E-4;
% P = D;          %Jacobi
% Q = -(L+U);     %Jacobi
% 
% for t = 1:100
%     d = Q*x + b;
%     for i=1:n
%         x(i) = d(i)/P(i,i);
%     end
%     if norm(b - A*x)<tol
%         return
%     end
% end
    

% %Gauss-Seidel
% tol = 1E-4;
% P = L+D;        %Gauss-Seidel 
% Q = -U;         %Gauss-Seidel
% 
% for t = 1:100
%     x = Q*x + b;
%     d = x;
%     for i = 1:n
%         x(i) = (x(i) - P(i,1:i-1)*x(1:i-1))/P(i,i);
%     end
%     if norm(b - A*x)<tol
%         return
%     end
% end

% %SOR and Gauss-Seidel CSC
% 
% tol =1E-4;
% w = 1;
% 
% %Matrices
% m = length(Ac)-1;
% nzp = m + (length(Av)-m)/2;
% if w==1
%     nzq = (length(Av)-m)/2;
% else
%     nzq = m + (length(Av)-m)/2;
% end
% 
% Qv = zeros(nzq,1);
% Qr = zeros(nzq,1);
% Qc = zeros(m+1,1);
% Qc(1) = 1;
% Pv = zeros(nzp,1);
% Pr = zeros(nzp,1);
% Pc = zeros(m+1,1);
% Pc(1) = 1;
%     
% qp=1;
% pp=1;
% if w==1
%     for j=1:m
%         for i=Ac(j):Ac(j+1)-1
%             if Ar(i) < j
%                 Qv(qp) = -Av(i);
%                 Qr(qp) = Ar(i);
%                 qp = qp + 1;
%             else
%                 Pv(pp) = Av(i);
%                 Pr(pp) = Ar(i);
%                 pp = pp + 1;
%             end
%         end
%         Qc(j+1) = qp;
%         Pc(j+1) = pp;
%     end
% else
%     for j=1:m
%         for i=Ac(j):Ac(j+1)-1
%             if Ar(i) < j
%                 Qv(qp) = -Av(i);
%                 Qr(qp) = Ar(i);
%                 qp = qp + 1;
%             elseif Ar(i) > j
%                 Pv(pp) = Av(i);
%                 Pr(pp) = Ar(i);
%                 pp = pp + 1;
%             else
%                 Qv(qp) = (1/w - 1)*Av(i);
%                 Qr(qp) = Ar(i);
%                 qp = qp + 1;
%                 Pv(pp) = (1/w)*Av(i);
%                 Pr(pp) = Ar(i);
%                 pp = pp + 1;
%             end
%         end
%         Qc(j+1) = qp;
%         Pc(j+1) = pp;
%     end
% end
%     
% % loop
% for t = 1:100
%     x = csc_matvec(Qv,Qr,Qc,x) + b;
%     d = x;
%     for j=1:m
%         for i = Pc(j):Pc(j+1)-Pc(1)
%             if j==Pr(i)
%                 x(j) = x(j)/Pv(i);
%             else
%                 x(Pr(i)) = x(Pr(i)) - Pv(i)*x(j);
%             end
%         end
%     end
%     if norm(b - csc_matvec(Av,Ar,Ac,x))<tol
%         return
%     end
% end

% %SOR
% % if w = 1, Gauss-Seidel
% tol = 1E-4;
% w = 1;                %SOR 1.09 to 1.16 with tol = 1E-16
% P = (1/w)*D + L;        %SOR
% Q = (1/w - 1)*D - U;    %SOR
% 
% for t = 1:100
%     x = Q*x + b;
%     d = x;
%     for i = 1:n
%         x(i) = (x(i) - P(i,1:i-1)*x(1:i-1))/P(i,i);
%     end
%     if norm(b - A*x)<tol
%         return
%     end
% end

% % CG CSC
% tol = 1E-4;
% ro = b - csc_matvec(Av,Ar,Ac,x);
% d = ro;
% for t = 1:100
%     Ad = csc_matvec(Av,Ar,Ac,d);
%     alf = (ro'*ro)/(d'*Ad);
%     x = x + alf*d;
%     if norm(b - csc_matvec(Av,Ar,Ac,x))<tol
%         return
%     end
%     r = ro - alf*Ad;
%     bet = (r'*r)/(ro'*ro);
%     d = r + bet*d;
%     ro = r;
% end

% %CG
% tol = 14E-4;
% ro = b - A*x;
% d = ro;
% for t = 1:100
%     alf = (ro'*ro)/(d'*A*d);
%     x = x + alf*d;
%     if norm(b - A*x)<tol
%         return
%     end
%     r = ro - alf*A*d;
%     bet = (r'*r)/(ro'*ro);
%     d = r + bet*d;
%     ro = r;
% end

% % Preconditioned CG CSC
% m = length(Ac)-1;
% 
% prec = 2;
% 
% if prec==3
%     m = length(Ac)-1;
%     nzp = m + (length(Av)-m)/2;
%     Pv = zeros(nzp,1);
%     Pr = zeros(nzp,1);
%     Pc = zeros(m+1,1);
%     Pc(1) = 1;
%     pp = 1;
%     for j=1:m
%         for i=Ac(j):Ac(j+1)-1
%             if Ar(i) >= j
%                 Pv(pp) = Av(i);
%                 Pr(pp) = Ar(i);
%                 pp = pp + 1;
%             end
%         end
%         Pc(j+1) = pp;
%     end
% else
%     P = zeros(m,1);
%     for j=1:m
%         for i=Ac(j):Ac(j+1)-1
%             if Ar(i) == j
%                 P(j) = Av(i);
%             end
%         end
%     end
% end
% 
% tol = 1E-4;
% ro = b - csc_matvec(Av,Ar,Ac,x);
% 
% if prec==3
%     zo = ro;
%     for j=1:m
%         for i = Pc(j):Pc(j+1)-Pc(1)
%             if j==Pr(i)
%                 zo(j) = zo(j)/Pv(i);
%             else
%                 zo(Pr(i)) = zo(Pr(i)) - Pv(i)*zo(j);
%             end
%         end
%     end
% else
%     z = zeros(n,1);
%     zo = ro./P;
% end
% 
% d = zo;
% for t = 1:100
%     Ad = csc_matvec(Av,Ar,Ac,d);
%     alf = (ro'*zo)/(d'*Ad);
%     x = x + alf*d;
%     if norm(b - csc_matvec(Av,Ar,Ac,x))<tol
%         return
%     end
%     r = ro - alf*Ad;
% 
%     if prec==3
%         z = r;
%         for j=1:m
%             for i = Pc(j):Pc(j+1)-Pc(1)
%                 if j==Pr(i)
%                     z(j) = z(j)/Pv(i);
%                 else
%                     z(Pr(i)) = z(Pr(i)) - Pv(i)*z(j);
%                 end
%             end
%         end
%     else
%         z = r./P;
%     end
%     
%     bet = (r'*z)/(ro'*zo);
%     d = z + bet*d;
%     ro = r;
%     zo = z;
% end

% %Preconditioned CG
% P = D;
% % w = 1;                %SOR 1.09 to 1.16 with tol = 1E-16
% % P = (1/w)*D + L;        %SOR
% tol = 14E-4;
% ro = b - A*x;
% 
% % z = zeros(n,1);
% % zo = zeros(n,1);
% % for i=1:n
% %     zo(i) = ro(i)/P(i,i);
% % end
% 
% zo = ro;
% for i = 1:n
%     zo(i) = (zo(i) - P(i,1:i-1)*zo(1:i-1))/P(i,i);
% end
% 
% d = zo;
% for t = 1:100
%     alf = (ro'*zo)/(d'*A*d);
%     x = x + alf*d;
%     if norm(b - A*x)<tol
%         return
%     end
%     r = ro - alf*A*d;
% 
% %     for i=1:n
% %         z(i) = r(i)/P(i,i);
% %     end
%     
%     z = r;
%     for i = 1:n
%         z(i) = (z(i) - P(i,1:i-1)*z(1:i-1))/P(i,i);
%     end
%     
%     bet = (r'*z)/(ro'*zo);
%     d = z + bet*d;
%     ro = r;
%     zo = z;
% end

% %matmul and transposition
% A = [0 0 2 0 -1 0;4 0 3 3 7 0;-2 0 0 0 0 -1; 0 1 0 1 0 0];
% B = [1 0 -1 0 0 5;0 0 0 0 -2 0;4 6 0 2 0 0;0 -1 1 0 0 0];
% [Av,Ac,Ar] = full2csr(A);
% [Bv,Bc,Br] = full2csr(B);
% [Bv,Bc,Br] = css_trans(Bv,Bc,Br);
% 
% [Cc,Cr] = csr_symmatmul(Ac,Ar,Bc,Br);
% [Cc,Cr] = css_symtrans(Cc,Cr);
% [Cc,Cr] = css_symtrans(Cc,Cr);
% [Cv,Cc,Cr] = csr_nummatmul(Av,Ac,Ar,Bv,Bc,Br,Cc,Cr);

% %Kronecker product
% A = [1 1 0; 0 1 0;0 0 2];
% B = [1 2; 3 4];
% C = kron(A,B);
% [Acv,Acr,Acc] = full2csc(A);
% [Bcv,Bcr,Bcc] = full2csc(B);
% [Arv,Arc,Arr] = full2csr(A);
% [Brv,Brc,Brr] = full2csr(B);
% [Ccv,Ccr,Ccc] = css_kron(Acv,Acr,Acc,Bcv,Bcr,Bcc);
% [Crv,Crc,Crr] = css_kron(Arv,Arc,Arr,Brv,Brc,Brr);
% % csc_kron works either for CSC and CSR

% % csr Jacobi
% A = [1 0 0; 0 2 0;0 0 3];
% B = [3 1; 1 4];
% C = kron(A,B);
% [Cv,Cc,Cr] = full2csr(C);
% [Pv,Qv,Qc,Qr] = csr_prejacobi(Cv,Cc,Cr);
% b = ones(6,1);
% % d = csr_matvec(Cv,Cc,Cr,b);
% niter = 20;
% tol = 1e-8;
% [xr,tr] = csr_jacobi(Cv,Cc,Cr,Pv,Qv,Qc,Qr,b,b,niter,tol);
% 
% [Ccv,Ccr,Ccc] = full2csc(C);
% [Pcv,Qcv,Qcr,Qcc] = csc_prejacobi(Ccv,Ccr,Ccc);
% % d = csc_matvec(Ccv,Ccr,Ccc,b)
% [xc,tc] = csc_jacobi(Ccv,Ccr,Ccc,Pcv,Qcv,Qcr,Qcc,b,b,niter,tol);
% 
% x = C\b;

% % csr SOR
% Aa = [1 0 0; 0 2 0;0 0 3];
% Ab = [3 1; 1 4];
% A = kron(Aa,Ab);
% A(3,1) = 1;
% A(3,2) = 1;
% A(4,1) = 1;
% [Av,Ac,Ar] = full2csr(A);
% [Acv,Acr,Acc] = full2csc(A);
% w = 1;
% [Pv,Pc,Pr,Qv,Qc,Qr] = csr_preSOR(Av,Ac,Ar,w);
% [Pcv,Pcr,Pcc,Qcv,Qcr,Qcc] = csc_preSOR(Acv,Acc,Acr,w);
% % b = ones(6,1);
% b = [1;2;3;4;5;6];
% niter = 20;
% tol = 1e-8;
% % x = A\b;
% [xr,tr] = csr_SOR(Av,Ar,Ac,Pv,Pr,Pc,Qv,Qr,Qc,b,b,niter,tol);
% [xc,tc] = csc_SOR(Acv,Acr,Acc,Pcv,Pcr,Pcc,Qcv,Qcr,Qcc,b,b,niter,tol);

% % %csr Conjugate Gradient
% Aa = [1 0 0; 0 2 0;0 0 3];
% Ab = [3 1; 1 4];
% A = kron(Aa,Ab);
% [Av,Ac,Ar] = full2csr(A);
% b = [1;2;3;4;5;6];
% x = A\b;
% niter = 20;
% tol = 1e-8;
% [xr,t] = csr_CG(Av,Ac,Ar,b,b,niter,tol);
% [xrp,tp] = csr_DPCG(Av,Ac,Ar,b,b,niter,tol);

% % %Projection Methods Experiments (Gauss-Seidel)
% Aa = [1 0 0; 0 2 0;0 0 3];
% Ab = [3 1; 1 4];
% A = kron(Aa,Ab);
% [Av,Ac,Ar] = full2csr(A);
% [Pv,Pc,Pr,Qv,Qc,Qr] = csr_preSOR(Av,Ac,Ar,1); %Gauss-Seidel
% b = [1;1;1;1;1;1];
% x0 = [0;1;0;0;0;0];
% niter = 20;
% tol = 1e-8;
% xm = A\b;
% 
% n = length(Ar)-1;
% ro = b - csr_matvec(Av,Ac,Ar,x0);
% 
% % loop
% % for t = 1:niter
%     xg = csr_matvec(Qv,Qc,Qr,x0) + b;
%     for i=1:n
%         j = Pr(i):Pr(i+1)-Pr(1);
%         if length(j)==1
%             xg(i) = xg(i)/Pv(j(length(j)));
%         else
%             xg(i) = (xg(i) - Pv(j(1:length(j)-1))*xg(Pc(j(1:length(j)-1))))/Pv(j(length(j)));
%         end
%     end
%     rnew = b - csr_matvec(Av,Ac,Ar,xg);
%     xg'*rnew
%     d = xg - x0;
%     rnew'*d
%     wo = ro - csr_matvec(Av,Ac,Ar,d);
%     d'*wo
%     x0 = xg;
% %     if norm(b - csr_matvec(Av,Ac,Ar,x))<tol
% %         return
% %     end
% % end

% %Projection Methods Experiments (This is Gauss-Seidel in Projection Form)
% h = [1 3 1;1 3 1;2 6 2;2 6 2;3 8 3;3 8 3];
% A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
% [Av,Ac,Ar] = full2csr(A);
% b = ones(length(A),1);
% x0 = ones(length(A),1);
% niter = 20;
% tol = 1e-8;
% xm = A\b;
% % err = zeros(1,100);
% normr = zeros(1,length(A));
% V = A*b;
% V = V./(norm(V));
% 
% for i = 1:length(A)
%     r = b - csr_matvec(Av,Ac,Ar,x0);
%     %space K
%     %
% %     V=zeros(length(A),i);
% %     for j=1:i
% %         V(j,j) = 1;
% %     end
% %     V(1) = 1;
%     %
% %     V(1 + mod(i,6),2) = 1;
%     %
% %     V = ones(6,1);
% %     V(1 + mod(i-1,6),1) = 0;
%     %
%     V = [r A*r A*A*r];    %krylov subspace
% 
% %     space L
%     W=A*V;    %Petrov-Galerkin
% %     W=V;        %Galerkin
%     
%     WAV = W'*A*V;
%     y = inv(WAV)*W'*r;
%     xg = x0 + V*y;
%     d = xg - x0;
%     w0 = r - A*d;
%     ort = (W'*w0)'
%     normr(i) = norm(b - csr_matvec(Av,Ac,Ar,xg));
%     x0 = xg;
% %     err(i) = norm(xm-x0);
% 
% %
% %     Va = A*V(:,end);
% %     Va = Va./(norm(Va));
% %     V = [V Va];
% % V = [V(:,end) Va];
% 
% end
% % plot(err)
% % ylim([0 3])
% %of course the selection of basis for K and L are way important for the
% %convergence!!!








