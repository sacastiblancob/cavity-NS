function [ABv,ABr,ABc] = csc_sum(Av,Ar,Ac,Bv,Br,Bc)
%Computes the sum of A+B with A and B stored in CSC
% Av and Bv --> values of A and B in CSC respectively
% Br and Br --> rows indices of A and B in CSC respectively
% Ac and Bc --> columns indices of A and B in CSC respectively
%
% Returns ABv, ABr and ABc with AB=A+B
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%      Pontificia Universidad Javeriana - Bogotá
%

ma = length(Ac)-1;
mb = length(Bc)-1;
if (ma~=mb)
    disp('Dimensions does not agree!!!, errors might be expected!!!')
    ABv = NaN;
    ABr = NaN;
    ABc = NaN;
    return
else
    m = ma;
end
%nza = Ac(m+1)-1;
%nzb = Bc(m+1)-1;

% if length(Av)==length(Bv) && length(Ar)==length(Br) && length(Ac)==length(Bc)
%     if min(Av==Bv)==1 && min(Ar==Br)==1 && min(Ac==Bc)==1
%         ABv = Av.*2;
%         ABr = Ar;
%         ABc = Ac;
%         return
%     elseif min(Av==-Bv)==1 && min(Ar==Br)==1 && min(Ac==Bc)==1
%         ABv = 0;
%         ABr = 0;
%         ABc = 0;
%         return
%     end
% end
    
ABc = ones(m+1,1);
nzab = 0;
for j=1:m
    ra = Ar(Ac(j):Ac(j+1)-Ac(1));
    rb = Br(Bc(j):Bc(j+1)-Bc(1));
    r = union(ra,rb);
    va = Av(Ac(j):Ac(j+1)-Ac(1));
    vb = Bv(Bc(j):Bc(j+1)-Bc(1));
    for i = r'
            if (ismember(i,ra) && ismember(i,rb))
                if (va(ra==i)+vb(rb==i)==0)
                    continue
                else
                    nzab = nzab + 1;
                    ABc(j+1:m+1) = ABc(j+1:m+1) + 1;
                end
            else
                nzab = nzab + 1;
                ABc(j+1:m+1) = ABc(j+1:m+1) + 1;
            end
    end
end
ABv = zeros(nzab,1);
ABr = zeros(nzab,1);
k=1;
for j=1:m
    ra = Ar(Ac(j):Ac(j+1)-Ac(1));
    rb = Br(Bc(j):Bc(j+1)-Bc(1));
    r = union(ra,rb);
    va = Av(Ac(j):Ac(j+1)-Ac(1));
    vb = Bv(Bc(j):Bc(j+1)-Bc(1));
    for i = r'
        if (ismember(i,ra) && ismember(i,rb))
            if (va(ra==i)+vb(rb==i)==0)
                continue
            else
                ABv(k) = va(ra==i)+vb(rb==i);
                ABr(k) = i;
                k=k+1;
            end
        else
            if ismember(i,ra)
                ABv(k) = va(ra==i);
                ABr(k) = i;
                k=k+1;
            else
                ABv(k) = vb(rb==i);
                ABr(k) = i;
                k=k+1;
            end
        end
    end
end

end