function [liam,ljbl] = cross_lagrange2(ep,X2,et,Y2,Nx,Ny)
%
liam = zeros(Nx+1,Ny+1);
ljbl = zeros(Ny+1,Nx+1);

for i=1:Ny+1
    for m=1:Nx+1
%         liam(m,i) = lagrange_point2(ep(m),X2(i,:));
        liam(m,i) = lagrange_point(ep(m),X2(i,:),i);
    end
end

for j=1:Nx+1
    for l=1:Ny+1
%         ljbl(l,j) = lagrange_point2(et(l),Y2(:,j));
        ljbl(l,j) = lagrange_point(et(l),Y2(:,j),j);
    end
end

end