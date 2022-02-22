function [LUv,LUr,LUc] = csc_packlu(Av,Ar,Ac)
%Function for compute the packed-LU decomposition for a matriz A stored in
%CSC (Compressed Sparse Column, beggining counting in 1 not in 0)
%Entries:
%       Av = Vector with values of A
%       Ar = Vector with rows indexes of A in CSC format
%       Ac = Vector with columns indexes of A in CSC format
%

m = length(Ac)-Ac(1);
nz = length(Av);

LUv = Av;
LUr = Ar;
LUc = Ac;

%k=1;
fc=0;       %value of the column for which one should divide every step of LU
nzp = 0;    %non-zeros padded, if one zero becomes a non-zero one ought insert that element
            %col = column identifier for internal loop for variable l
            %fc = value that is used for operate the part that should be operated
            %  for example, if I'm zeroing out the A(2,1) entry, fc gonna
            %  be the entries A(1,*), the value of * depends where the
            %  loops are. If A(2,1) is non-zero, and supose that A(2,2) is
            %  also a non-zero, then fc=A(1,2), and then is used for
            %  compute LU(2,2)  = A(2,2) - (A(2,1)/A(1,1))*A(1,2)
            %                   = A(2,2) - (A(2,1)/d)*A(1,2)
            %Then d = value in the diagonal that correspond with where loop
            %that's it: the value in the diagonal that is used for
            %transform in zero the low part of A

%tic
%tol = 1E-10;
for j=1:m
    for i = LUc(j):LUc(j+1)-LUc(1)
        if j==LUr(i)    %this if proof if the value is on the diagonal
            d=LUv(i);
        end
        if j<LUr(i)             %this if proofs if the value is under the diagonal
            LUv(i) = LUv(i)/d;  %here is computed the new value if is under the diagonal, because need to become 0 in U and LUv(i)/d in L
            col = j;            %this is the actualization of col, which begins with the value of j
            for l = (i+1):nz    %this for walk the rest of the vectors, that's it: l always is bigger than i
                if (LUr(l+nzp) > LUr(i) && fc~=0) || (LUr(l+nzp)<LUr(l+nzp-1) && fc~=0)     %this if proof if a zero becomes a non-zero, it may be the case
                    val = LUv(l+nzp) - LUv(i)*fc;
%                     if abs(val) < tol
%                         break
%                     end
                    LUv = [LUv(1:l+nzp-1);0;LUv(l+nzp:nz+nzp)];                                 %here a zero is padded in vector of values v
                    LUv(l+nzp) = val;                                    %here the padded zero becomes the non-zero
                    
                    LUr = [LUr(1:l+nzp-1);LUr(i);LUr(l+nzp:nz+nzp)];                            %here the new non-zero row index is padded into the vector of row indexes r
                    
                    LUc = [LUc(1:col);(LUc(col+1:end)+1)];                                  %here the column indexes bigger than c(col)increase by one (the effect of padding a nes non-zero)
                    
                    nzp = nzp+1;                                                            %here nzp increase 1, because a new non-zero was padded
                    fc=0;                                                                   %column value for zeroed the Low part of U is zeroed
                end
                if j==LUr(l+nzp)                                                            %this if proofs if l arrives to the diagonal (which is the value of fc)
                    fc = LUv(l+nzp);                                                        %takes the column value for zeroed the Low part of U
                end
                if LUr(l+nzp) == LUr(i)                                                     %this if proofs if l arrives to a value in the same row than the value that ones is zeroed out in U
                    LUv(l+nzp) = LUv(l+nzp) - LUv(i)*fc;                                    %here is computed the typical LU operation
                    fc = 0;                                                                 %column value for zeroed the Low part of U is zeroed
                end
                
                if LUr(l+nzp)<=LUr(l+nzp-1)                                                  %this if proofs if one change of column
                    col=col+1;                                                              %as one has change of column, here col increase by one
                end
                if l==nz                                                                    %this if proofs if one arrives to the end of the vectors v and r
                    nz=nz + nzp;                                                            %the number of non-zeros(nz) is updated with the nzp value
                    nzp=0;                                                                  %the non-zero padded value is zeroed
                end
            end           
        end
    end
end
%toc

% %proof for the case were vector v keep zeros
% tol = 1E-16;
% col = 1;
% zf = 0;
% for i=1:nz
%     if (i-zf)>1
%         if LUr(i-zf) <= LUr(i-zf-1)
%             col = col + 1;
%         end
%     end
%     if abs(LUv(i-zf)) <= tol
%         LUv = [LUv(1:i-zf-1);LUv(i-zf+1:nz)];
%         LUr = [LUr(1:i-zf-1);LUr(i-zf+1:nz)];
%         LUc = [LUc(1:col);(LUc(col+1:end)-1)];
%         nz = nz-1;
%         zf = zf + 1;
%     end
% end

end