function [v,beta] = householderv(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the vector v for constructing Householder
% Projector related to vector x, i.e
%           P = I - beta*v*v'
% Then,          Px = ||x||*e_1
%
% Algorithm 5.1.1 (Householder Vector) Matrix Computations (Golub & Van Loan)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(x);
x = x/norm(x);      %to avoid overflow
sig = x(2:n)'*x(2:n);
v = [1;x(2:n)];
if n==1
    beta = 2;
    v = 1;
    return
end
if sig==0
    beta=0;
else
    mu = sqrt(x(1)^2+sig);
    if x(1) <= 0
        v(1) = x(1) - mu;
    else
        v(1) = -sig/(x(1) + mu);
    end
    beta = (2*v(1)^2)/(sig+v(1)^2);
    v = v/v(1);
end

end