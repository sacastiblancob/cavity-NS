function csc_spy(Bv,Br,Bc)
%
% Spy plot of a CSC stored matrix
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%      Pontificia Universidad Javeriana - Bogotá
%

nz = length(Bv);
cols = zeros(nz,1);
m = length(Bc)-1;
n = max(Br);
k=1;
for j=1:m
    for i=Bc(j):Bc(j+1)-1
        cols(k) = j;
        k = k + 1;
    end
end

figure
h1 = axes;
scatter(cols,Br)

ylim([0.5 n+0.5])
set(h1,'Ydir','reverse')
xlabel(['nz = ' num2str(nz)])
axis equal
xlim([0.5 m+0.5])

end
