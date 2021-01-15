function [U,V] = deriv2D_up2(Uo,Vo,X,Y,dt)
    dx = abs(X(1,2)-X(1,1));
    dy = abs(Y(2,1)-Y(1,1));
    sizex = size(X);
    U = zeros(size(Uo));
    V = zeros(size(Vo));
    m = sizex(1);
    n = sizex(2);
    vecpos = 1:n*m';
    Mpos = reshape(vecpos,n,m);
%     bound = [Mpos(1,:), Mpos(n,:), Mpos(2:n-1,1)',Mpos(2:n-1,m)'];
%     bound = sort(bound);
    boundint = [Mpos(2,2:m-1), Mpos(n-1,2:m-1), Mpos(3:n-2,2)',Mpos(3:n-2,m-1)'];
    boundint = sort(boundint);
    for i = boundint
        if Uo(i) >= 0 && Vo(i) >= 0
            U(i) = Uo(i) - (dt/dx)*Uo(i)*(Uo(i) - Uo(i-n)) - (dt/dy)*Vo(i)*(Uo(i) - Uo(i-1));
            V(i) = Vo(i) - (dt/dx)*Uo(i)*(Vo(i) - Vo(i-n)) - (dt/dy)*Vo(i)*(Vo(i) - Vo(i-1));
        elseif Uo(i) > 0 && Vo(i) < 0
            U(i) = Uo(i) - (dt/dx)*Uo(i)*(Uo(i) - Uo(i-n)) - (dt/dy)*Vo(i)*(Uo(i+1) - Uo(i));
            V(i) = Vo(i) - (dt/dx)*Uo(i)*(Vo(i) - Vo(i-n)) - (dt/dy)*Vo(i)*(Vo(i+1) - Vo(i));
        elseif Uo(i) < 0 && Vo(i) > 0
            U(i) = Uo(i) - (dt/dx)*Uo(i)*(Uo(i+n) - Uo(i)) - (dt/dy)*Vo(i)*(Uo(i) - Uo(i-1));
            V(i) = Vo(i) - (dt/dx)*Uo(i)*(Vo(i+n) - Vo(i)) - (dt/dy)*Vo(i)*(Vo(i) - Vo(i-1));
        else
            U(i) = Uo(i) - (dt/dx)*Uo(i)*(Uo(i+n) - Uo(i)) - (dt/dy)*Vo(i)*(Uo(i+1) - Uo(i));
            V(i) = Vo(i) - (dt/dx)*Uo(i)*(Vo(i+n) - Vo(i)) - (dt/dy)*Vo(i)*(Vo(i+1) - Vo(i));
        end
    end
    for i = 3:m-2
        for j = 3:n-2
            if Uo(i,j) >= 0 && Vo(i,j) > 0
                U(i,j) = Uo(i,j) - (dt/(2*dx))*Uo(i,j)*(3*Uo(i,j) - 4*Uo(i,j-1) + Uo(i,j-2)) - (dt/(2*dy))*Vo(i,j)*(3*Uo(i,j) - 4*Uo(i-1,j) + Uo(i-2,j));
                V(i,j) = Vo(i,j) - (dt/(2*dx))*Uo(i,j)*(3*Vo(i,j) - 4*Vo(i,j-1) + Vo(i,j-2)) - (dt/(2*dy))*Vo(i,j)*(3*Vo(i,j) - 4*Vo(i-1,j) + Vo(i-2,j));
            elseif Uo(i,j) >= 0 && Vo(i,j) <= 0
                U(i,j) = Uo(i,j) - (dt/(2*dx))*Uo(i,j)*(3*Uo(i,j) - 4*Uo(i,j-1) + Uo(i,j-2)) - (dt/(2*dy))*Vo(i,j)*(-Uo(i+2,j) + 4*Uo(i+1,j) - 3*Uo(i,j));
                V(i,j) = Vo(i,j) - (dt/(2*dx))*Uo(i,j)*(3*Vo(i,j) - 4*Vo(i,j-1) + Vo(i,j-2)) - (dt/(2*dy))*Vo(i,j)*(-Vo(i+2,j) + 4*Vo(i+1,j) - 3*Vo(i,j));
            elseif Uo(i,j) < 0 && Vo(i,j) > 0
                U(i,j) = Uo(i,j) - (dt/(2*dx))*Uo(i,j)*(-Uo(i,j+2) + 4*Uo(i,j+1) - 3*Uo(i,j)) - (dt/(2*dy))*Vo(i,j)*(3*Uo(i,j) - 4*Uo(i-1,j) + Uo(i-2,j));
                V(i,j) = Vo(i,j) - (dt/(2*dx))*Uo(i,j)*(-Vo(i,j+2) + 4*Vo(i,j+1) - 3*Vo(i,j)) - (dt/(2*dy))*Vo(i,j)*(3*Vo(i,j) - 4*Vo(i-1,j) + Vo(i-2,j));
            elseif Uo(i,j) < 0 && Vo(i,j) <= 0
                U(i,j) = Uo(i,j) - (dt/(2*dx))*Uo(i,j)*(-Uo(i,j+2) + 4*Uo(i,j+1) - 3*Uo(i,j)) - (dt/(2*dy))*Vo(i,j)*(-Uo(i+2,j) + 4*Uo(i+1,j) - 3*Uo(i,j));
                V(i,j) = Vo(i,j) - (dt/(2*dx))*Uo(i,j)*(-Vo(i,j+2) + 4*Vo(i,j+1) - 3*Vo(i,j)) - (dt/(2*dy))*Vo(i,j)*(-Vo(i+2,j) + 4*Vo(i+1,j) - 3*Vo(i,j));
            end
        end
    end
end