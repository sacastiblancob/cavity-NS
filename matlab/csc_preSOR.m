function [Pv,Pr,Pc,Qv,Qr,Qc] = csc_preSOR(Av,Ar,Ac,w)
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
%     Av, Ar, Ac : Components of matrix A in CSC stotage
%     w : SOR Over-Relaxation coefficient
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%      Pontificia Universidad Javeriana - Bogotá
%

m = length(Ac)-1;
nzp = 0;
nzq = 0;

for j=1:m
    for i=Ac(j):Ac(j+1)-1
        if Ar(i)<j
            nzq = nzq + 1;
        elseif Ar(i)>j
            nzp = nzp + 1;
        elseif Ar(i)==j
            if w==1
                nzp = nzp+1;
            else
                nzp = nzp+1;
                nzq = nzq+1;
            end
        end
    end
end

Qv = zeros(nzq,1);
Qr = zeros(nzq,1);
Qc = zeros(m+1,1);
Qc(1) = 1;
Pv = zeros(nzp,1);
Pr = zeros(nzp,1);
Pc = zeros(m+1,1);
Pc(1) = 1;
    
qp=1;
pp=1;
if w==1
    for j=1:m
        for i=Ac(j):Ac(j+1)-1
            if Ar(i) < j
                Qv(qp) = -Av(i);
                Qr(qp) = Ar(i);
                qp = qp + 1;
            else
                Pv(pp) = Av(i);
                Pr(pp) = Ar(i);
                pp = pp + 1;
            end
        end
        Qc(j+1) = qp;
        Pc(j+1) = pp;
    end
else
    for j=1:m
        for i=Ac(j):Ac(j+1)-1
            if Ar(i) < j
                Qv(qp) = -Av(i);
                Qr(qp) = Ar(i);
                qp = qp + 1;
            elseif Ar(i) > j
                Pv(pp) = Av(i);
                Pr(pp) = Ar(i);
                pp = pp + 1;
            else
                Qv(qp) = (1/w - 1)*Av(i);
                Qr(qp) = Ar(i);
                qp = qp + 1;
                Pv(pp) = (1/w)*Av(i);
                Pr(pp) = Ar(i);
                pp = pp + 1;
            end
        end
        Qc(j+1) = qp;
        Pc(j+1) = pp;
    end
end

end