function [H,Q] = hessenberg(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: Compute Hessenberg similar matrix to A
%
% Entries: A
%
% Output: H,Q
% Such as:
%
%         A = Q'*H*Q
%
%         H ---> Hessenberg matrix
%         Q ---> Ortogonal matrix
%
% If A is symmetric, H is tridiagonal (-1,0,1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(A);
Q = eye(n);
H = A;

for i=2:n
    u = H(i:n,i-1);
    if norm(u)~=0
        if u(1) ~= 0
           u(1) = u(1) + (u(1)/abs(u(1)))*norm(u);
        else
           u(1) = u(1) + norm(u);
        end
        u = u./norm(u);
    end
    
    P = eye(n);
    P(i:n,i:n) = eye(n-i+1) - 2.*u*u';
    Q = P*Q;
    H = P*H*P;
    i
end

end
