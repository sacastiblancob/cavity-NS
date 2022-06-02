function [FI] = lagrange_2D_interpolation(ep,et,X2,Y2,F2,Nx,Ny)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function interpolates the values of F2(X2,Y2) into the coordinates
% ep, et ,i.t., FI(ep,et). Tipically (ep,et) are the coordinates of
% quadrature nodes in the master element.
%
%
% Based on: DG-FTLE: Lagrangian coherent structures with high-order
%           discontinuous-Galerkin methods
%               Daniel A. Nelson, Gustaaf B. Jacobs, 2015
% &
% David Kopriva. Spectral Element Methods
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ILA = zeros(Ny+1,(Nx+1)*(Ny+1));
ILB = zeros((Nx+1)*(Ny+1),Nx+1);

for j=1:(Nx+1)
    ILA(:,(j-1)*(Nx+1)+1:j*(Nx+1)) = cross_lagrange_ep(ep(j),X2,Nx,Ny);
end
for i=1:(Ny+1)
    ILB((i-1)*(Ny+1)+1:i*(Ny+1),:) = cross_lagrange_et(et(i),Y2,Nx,Ny);
end

% psip = reshape(F2,1,[]);
% FI = zeros(Ny+1,Nx+1);
% for j=1:Nx+1
%     for i=1:Ny+1
%         ll = reshape(ILA(:,(j-1)*(Nx+1)+1:j*(Nx+1)).*ILB((i-1)*(Ny+1)+1:i*(Ny+1),:),1,[])';
% %         fi = psip*ll;
%         FI(i,j) = psip*ll;
%     end
% end

psip = reshape(F2,1,[]);
FI = F2;
for j=2:Nx
    for i=2:Ny
        ll = reshape(ILA(:,(j-1)*(Nx+1)+1:j*(Nx+1)).*ILB((i-1)*(Ny+1)+1:i*(Ny+1),:),1,[])';
%         fi = psip*ll;
        FI(i,j) = psip*ll;
    end
end

end

