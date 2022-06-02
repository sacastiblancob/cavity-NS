function D = dmatrix(xi,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes derivation matrix with Galerkin approach
% - Lagrange polynomials
%

ddiag = zeros(N+1,1);
for j=1:N+1
    aux = (xi(j) - xi);
    ddiag(j) = sum(1./aux(aux~=0));
end

D = diag(ddiag);

for j = 1:N+1
    for i = j+1:N+1
        num = (xi(i) - xi);
        den = (xi(j) - xi);
        kprod = num(den~=0)./den(den~=0);
        kprod = (1/(xi(j) - xi(i)))*prod(kprod(kprod~=0));
        D(i,j) = kprod;
        D(N+2-i,N+2-j) = -kprod;
    end
end

end