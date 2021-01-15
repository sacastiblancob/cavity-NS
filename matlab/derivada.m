function [DU] = derivada(P,dx,dy)

sizep = size(P);
n = sizep(1);
m = sizep(2);
DPx = sparse(size(P));
DPy = sparse(size(P));
vecpos = 1:n*m';
Mpos = reshape(vecpos,n,m);
boundupd = Mpos(1,2:m-1);
bounddownd = Mpos(n,2:m-1);
boundleftd = Mpos(2:n-1,1)';
boundrightd = Mpos(2:n-1,m)';
boundintd = [Mpos(2,2:m-1), Mpos(n-1,2:m-1), Mpos(3:n-2,2)',Mpos(3:n-2,m-1)'];
boundintd = sort(boundintd);
for i = 3:m-2
    for j = 3:n-2
        DPx(i,j) = (1/(12*dx))*(U(i,j-2) - 8*U(i,j-1) + 8*U(i,j+1) - U(i,j+2));
        DPy(i,j) = (1/(12*dy))*(V(i-2,j) - 8*V(i-1,j) + 8*V(i+1,j) - V(i+2,j));
    end
end
end