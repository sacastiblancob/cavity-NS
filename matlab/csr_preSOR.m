function [Pv,Pc,Pr,Qv,Qc,Qr] = csr_preSOR(Av,Ac,Ar,w)
%
% This function takes the matrix A in CSC storage, and prepares the
% matrices for solve the system Ax=b with SOR iterative method
%
% P and Q:
% Classic iterative methods use the next recurrence equation for solve the 
% system:
% 
%               P * x(t+1) = Q * x(t) + b
%
%       For SOR:
%            P = ((1/w)*D + L)
%            Q = ((1/w - 1)*D - U); where: L = lower part of A
%                                          U = upper part of A
%                                          D = diagonal of A
%
% For the particular case when w = 1, we are in Gauss-Seidel method
%
% Entries:
%     Av, Ar, Ac : Components of matrix A in CSR stotage
%     w : SOR Over-Relaxation coefficient
%
%%%%%%%

n = length(Ar)-1;
nzp = n;
if w==1
    nzq = 0;
else
    nzq = n;
end

for j=1:n
    for i=Ar(j):Ar(j+1)-1
        if Ac(i)<j
            nzq = nzq + 1;
        elseif Ac(i)>j
            nzp = nzp + 1;
        end
    end
end

Qv = zeros(1,nzq);
Qc = zeros(1,nzq);
Qr = zeros(1,n+1);
Qr(1) = 1;
Pv = zeros(1,nzp);
Pc = zeros(1,nzp);
Pr = zeros(1,n+1);
Pr(1) = 1;
    
qp=1;
pp=1;
if w==1
    for j=1:n
        for i=Ar(j):Ar(j+1)-1
            if Ac(i) > j
                Qv(qp) = -Av(i);
                Qc(qp) = Ac(i);
                qp = qp + 1;
            else
                Pv(pp) = Av(i);
                Pc(pp) = Ac(i);
                pp = pp + 1;
            end
        end
        Qr(j+1) = qp;
        Pr(j+1) = pp;
    end
else
    for j=1:n
        for i=Ar(j):Ar(j+1)-1
            if Ac(i) > j
                Qv(qp) = -Av(i);
                Qc(qp) = Ac(i);
                qp = qp + 1;
            elseif Ac(i) < j
                Pv(pp) = Av(i);
                Pc(pp) = Ac(i);
                pp = pp + 1;
            else
                Qv(qp) = (1/w - 1)*Av(i);
                Qc(qp) = Ac(i);
                qp = qp + 1;
                Pv(pp) = (1/w)*Av(i);
                Pc(pp) = Ac(i);
                pp = pp + 1;
            end
        end
        Qr(j+1) = qp;
        Pr(j+1) = pp;
    end
end

end