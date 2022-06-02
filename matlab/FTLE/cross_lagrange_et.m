function [lmbj] = cross_lagrange_et(et1,ett1,Nx,Ny)
%

lmbj = zeros(Ny+1,Nx+1);
for j=1:Nx+1
    for i=1:Ny+1
        lmbj(i,j) = lagrange_point(et1,ett1(:,j),i);
    end
end

end