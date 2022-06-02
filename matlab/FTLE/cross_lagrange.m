function [llai,lmbj] = cross_lagrange(ep1,ept1,et1,ett1,Nx,Ny)
%

llai = zeros(Ny+1,Nx+1);
for j=1:Nx+1
    for i=1:Ny+1
        llai(i,j) = lagrange_point(ep1,ept1(i,:),j);
    end
end

lmbj = zeros(Ny+1,Nx+1);
for j=1:Nx+1
    for i=1:Ny+1
        lmbj(i,j) = lagrange_point(et1,ett1(:,j),i);
    end
end

% llai = zeros(Nx+1,Ny+1);
% for j=1:Ny+1
%     for i=1:Nx+1
%         al = ep(i);
%         num = al - ept1(j,:);
%         den = ept1(j,j) - ept1(j,:);
%         llai(i,j) = prod(num(den~=0)./den(den~=0));
%     end
% end
% 
% lmbj = zeros(Ny+1,Nx+1);
% for j=1:Nx+1
%     for i=1:Ny+1
%         be = et(i);
%         num = be - ett1(:,j);
%         den = ept1(j,j) - ett1(:,j);
%         lmbj(i,j) = prod(num(den~=0)./den(den~=0));
%     end
% end

end