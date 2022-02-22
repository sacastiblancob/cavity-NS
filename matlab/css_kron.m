function [Cv,Cr,Cc] = css_kron(Av,Ar,Ac,Bv,Br,Bc)
%
% This function computes the kronecker product between A and B with A and B
% stored in CSC
%
% Example
%          A               B
%      |1  1  0|        |1 1 0|
%      |2  2  2|        |0 1 0|
%      |0  3  3|        |0 1 1|
%
% Output
%                      
%                     |1 1   1 1        |
%                     |  1     1        |
%          C          |  1 1   1 1      |
%      |1B  1B  0 |   |2 2   2 2   2 2  |
%      |2B  2B  2B| = |  2     2     2  |
%      |0   3B  3B|   |  2 2   2 2   2 2|
%                     |      3 3   3 3  |
%                     |        3     3  |
%                     |        3 3   3 3|
%                     
% The output, as the entries, is stored in CSC 
%
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%      Pontificia Universidad Javeriana - Bogotá
%

ma = length(Ac)-1;
mb = length(Bc)-1;
nza = Ac(ma+1)-1;
nzb = Bc(mb+1)-1;
nb = max(Br);
mc = ma*mb;
nzc = nza*nzb;
Cc = zeros(mc+1,1);
Cv = zeros(nzc,1);
Cr = zeros(nzc,1);
Cc(1) = 1;
aj=1;
bj=1;
p=1;
for j=1:mc
    af = (Ac(aj+1) - Ac(aj));
    bf = (Bc(bj+1) - Bc(bj));
    Cc(j+1) = Cc(j) + bf*af;
    r = nb;
    for ia = Ac(aj):Ac(aj+1)-1
        for ib = Bc(bj):Bc(bj+1)-1
            Cr(p) = r*(Ar(ia)-1)+Br(ib);
            Cv(p) = Av(ia)*Bv(ib);
            p = p+1;
        end
    end
    bj=bj+1;
    if bj==mb+1
        aj=aj+1;
        bj=1;
    end
end

end