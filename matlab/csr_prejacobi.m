function [Pv,Qv,Qc,Qr] = csr_prejacobi(Av,Ac,Ar)
%
% This function takes the matrix A in CSR storage, and prepares the
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
% Entries are the values(Av), colums(Ac) and rows-acum-index(Ar) in 
% CSR storage of A
%
%%%%

n = length(Ar)-1;
nzd = length(Av)-n;
Qv = zeros(1,nzd);
Qc = zeros(1,nzd);
Qr = zeros(1,n+1);
Qr(1) = 1;
Pv = zeros(1,n);
qp=1;
pp=1;
for i=1:n
    for j=Ar(i):Ar(i+1)-1
        if Ac(j) ~= i
            Qv(qp) = -Av(j);
            Qc(qp) = Ac(j);
            qp = qp + 1;
        else
            Pv(pp) = Av(j);
            pp = pp + 1;
        end
    end
    Qr(i+1) = qp;
end

end