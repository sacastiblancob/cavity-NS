function [Cv,Cr,Cc] = csc_matmul(Av,Ar,Ac,Bv,Br,Bc,epsi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the matrix-matrix multiplication with entry
% matrices stored in CSC
%
%          A*B = C
%
% epsi is the tolerance for zeros, if some value in the multiplication
% process becomes less than epsi, it will be taken as a computational zero
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dimensions of both entries matrices
ma = length(Ac)-1;
na = max(Ar);

mb = length(Bc)-1;
nb = max(Br);

% Proof about dimensions
if ma~=nb
    disp('Dimensions does not agree!!!, errors might be expected!!!')
    Cv = NaN;
    Cr = NaN;
    Cc = NaN;
    return
end

% Transposing A for better acces to information
[Tv,Tr,Tc] = csc_trans(Av,Ar,Ac);

mc = mb;
nzc = 0;
for bj=1:mb
    for aj=1:na
        ra = Tr(Tc(aj):Tc(aj+1)-Tc(1));
        rb = Br(Bc(bj):Bc(bj+1)-Bc(1));
        r = intersect(ra,rb);
        va = Tv(Tc(aj):Tc(aj+1)-Tc(1));
        vb = Bv(Bc(bj):Bc(bj+1)-Bc(1));
        v = 0;
        for i=1:length(r)
            v = v + va(ra==r(i))*vb(rb==r(i));
        end
        if abs(v) > epsi
            nzc = nzc + 1;
        end
    end
end

% Allocating memory for C
if nzc==0
    Cv = 0;
    Cr = 0;
    Cc = 0;
    return
else
    Cv = zeros(nzc,1);
    Cr = zeros(nzc,1);
    Cc = ones(mc+1,1);
end

% Filling C
k=1;
for bj=1:mb
    for aj=1:na
        ra = Tr(Tc(aj):Tc(aj+1)-Tc(1));
        rb = Br(Bc(bj):Bc(bj+1)-Bc(1));
        r = intersect(ra,rb);
        va = Tv(Tc(aj):Tc(aj+1)-Tc(1));
        vb = Bv(Bc(bj):Bc(bj+1)-Bc(1));
        v = 0;
        for i=1:length(r)
            v = v + va(ra==r(i))*vb(rb==r(i));
        end
        if abs(v) > epsi
            Cv(k) = v;
            Cr(k) = aj;
            k=k+1;
            Cc(bj+1:mc+1) = Cc(bj+1:mc+1) + 1;
        end
    end
end

end