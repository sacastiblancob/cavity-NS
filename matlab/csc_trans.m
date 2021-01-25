function [Tv,Tr,Tc] = csc_trans(Av,Ar,Ac)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes the matrix A in CSC storage, and computes its
% transpose such as:
%
%      T = A'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of columns and non-zeros of A
m = length(Ac)-1;
nz = Ac(m+1)-1;

%Number of columns of T
mt = max(Ar);

%Allocating Tv, Tr and Tc
Tv = zeros(nz,1);
Tr = zeros(nz,1);
Tc = zeros(mt+1,1);

% Making cols array
cols = zeros(nz,1);
for j=1:m
    for i=Ac(j):Ac(j+1)-1
        cols(i) = j;
    end
end

% Filling T arrays
Tc(1) = 1;
k=1;
for j=1:mt
    for i=1:nz    
        if Ar(i)==j
            Tr(k) = cols(i);
            Tv(k) = Av(i);
            Tc(j+1) = Tc(j+1)+1;
            k = k+1;
        end
    end
end

% Recursive computing over Tc
for j=2:mt+1
    Tc(j) = Tc(j) + Tc(j-1);
end

end