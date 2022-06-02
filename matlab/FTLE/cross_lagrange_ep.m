function [llai] = cross_lagrange_ep(ep1,ept1,Nx,Ny)
%

llai = zeros(Ny+1,Nx+1);
for j=1:Nx+1
    for i=1:Ny+1
        llai(i,j) = lagrange_point(ep1,ept1(i,:),j);
    end
end

end