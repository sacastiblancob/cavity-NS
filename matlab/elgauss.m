function x = elgauss(A,b)
%This function solves the system Ax=b by Gaussian Elimination

dim = size(A);
n = dim(1);
m = dim(2);
x = zeros(m,1);

%Computing the augmented matrix Ab
Ab = [A b];

%Computing the gaussian elimination over the augmented matrix Ab
for j = 1:n-1
    % Subtract multiples of row j to zero out Ab(i,j) for i > j
    Ab(j+1:n,:) = Ab(j+1:n,:) - Ab(j+1:n,j)/Ab(j,j)*Ab(j,:);
    
%     %null
%     Ab2(Ab2<tol) = 0;
%     U2 = Ab2(:,1:n);
%     spy(U2);
%     drawnow
%     pause(0.5)
end

% U = Ab(:,1:m);
% y = Ab(:,m+1);

%Performing the backward sustitution
for i = m:-1:1
    x(i) = (Ab(i,m+1) - Ab(i,i+1:m)*x(i+1:m))/Ab(i,i);
end

% for i = m:-1:1
%     if i==m
%         x(i) = 0;
%     elseif i==(m-1)
%         x(i) = 1;
%     else
%         x(i) = (Ab(i,m+1) - Ab(i,i+1:m)*x(i+1:m))/Ab(i,i);
%     end
% end

end