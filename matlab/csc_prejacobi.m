function [Pv,Qv,Qr,Qc] = csc_prejacobi(Av,Ar,Ac)
%
% This function takes the matrix A in CSC storage, and prepares the
% matrices for solve the system Ax=b with Jacobi iterative method
%
% P and Q:
% Classic iterative methods use the next recurrence equation for solve the 
% system:
% 
%               P * x(t+1) = Q * x(t) + b
%
%       For jacobi:
%            P = D = diagonal of matrix A
%            Q = -(L + U); where: L = lower part of A
%                                 U = upper part of A
%
% Entries are the values(Av), rows(Ar) and columns(Ac) in CSC storage of A
%
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%      Pontificia Universidad Javeriana - Bogotá
%

m = length(Ac)-1;
nzd = length(Av)-m;
Qv = zeros(nzd,1);
Qr = zeros(nzd,1);
Qc = zeros(m+1,1);
Qc(1) = 1;
Pv = zeros(m,1);
qp=1;
pp=1;
for j=1:m
    for i=Ac(j):Ac(j+1)-1
        if Ar(i) ~= j
            Qv(qp) = -Av(i);
            Qr(qp) = Ar(i);
            qp = qp + 1;
        else
            Pv(pp) = Av(i);
            pp = pp + 1;
        end
    end
    Qc(j+1) = qp;
end

end