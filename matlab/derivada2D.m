%Funcion para calcular las derivadas en dos dimensiones de un campo
%vectorial Uo y un campo vectorial Vo. El primero relacionado con la
%velocidad en X y el segundo con la velocidad en Y.
%Notacion:
%               [DUx,DVy,DUy,DVx] = derivada2D(Uo,Vo,dx,dy)
%Donde:               
%   Entrada:
%   Uo = campo de velocidades en X, matriz
%   Vo = campo de velocidades en Y, matriz
%   dx = diferencial numerico en direccion x, numero
%   dy = diferencial numerico en direccion y, numero
%   Salida:
%   DUx = matriz con las derivadas de Uo en direccion x
%   DVy = matriz con las derivadas de Vo en direccion y
%   DUy = matriz con las derivadas de Uo en direccion y
%   DVx = matriz con las derivadas de Vo en direccion x
%   
%   Nota: Se usa upwind de tres puntos para los nodos internos, upwind de
%   dos puntos para los nodos de la frontera y de los nodos vecinos a la
%   frontera. El metodo asegura solucion de adveccion si y solo si el
%   parametro de convergencia CFL <= 0.25

function [DUx,DVy,DUy,DVx] = derivada2D(Uo,Vo,dx,dy)
    sizex = size(Uo);
    DUx = zeros(size(Uo));
    DUy = zeros(size(Uo));
    DVx = zeros(size(Vo));
    DVy = zeros(size(Vo));
    m = sizex(1);
    n = sizex(2);
    vecpos = 1:n*m';
    Mpos = reshape(vecpos,n,m);
    boundupd = Mpos(1,2:m-1);
    bounddownd = Mpos(n,2:m-1);
    boundleftd = Mpos(2:n-1,1)';
    boundrightd = Mpos(2:n-1,m)';
    boundintd = [Mpos(2,2:m-1), Mpos(n-1,2:m-1), Mpos(3:n-2,2)',Mpos(3:n-2,m-1)'];
    boundintd = sort(boundintd);
    DUx(Mpos(1,1)) = (1/dx)*(Uo(Mpos(1,1)+n) - Uo(Mpos(1,1))); DUy(Mpos(1,1)) = (1/dy)*(Uo(Mpos(1,1)+1) - Uo(Mpos(1,1)));
    DVx(Mpos(1,1)) = (1/dx)*(Vo(Mpos(1,1)+n) - Vo(Mpos(1,1))); DVy(Mpos(1,1)) = (1/dy)*(Vo(Mpos(1,1)+1) - Vo(Mpos(1,1)));
    DUx(Mpos(1,n)) = (1/dx)*(Uo(Mpos(1,n)) - Uo(Mpos(1,n)-n)); DUy(Mpos(1,n)) = (1/dy)*(Uo(Mpos(1,n)+1) - Uo(Mpos(1,n)));
    DVx(Mpos(1,n)) = (1/dx)*(Vo(Mpos(1,n)) - Vo(Mpos(1,n)-n)); DVy(Mpos(1,n)) = (1/dy)*(Vo(Mpos(1,n)+1) - Vo(Mpos(1,n)));
    DUx(Mpos(m,1)) = (1/dx)*(Uo(Mpos(m,1)+n) - Uo(Mpos(m,1))); DUy(Mpos(m,1)) = (1/dy)*(Uo(Mpos(m,1)) - Uo(Mpos(m,1)-1));
    DVx(Mpos(m,1)) = (1/dx)*(Vo(Mpos(m,1)+n) - Vo(Mpos(m,1))); DVy(Mpos(m,1)) = (1/dy)*(Vo(Mpos(m,1)) - Vo(Mpos(m,1)-1));
    DUx(Mpos(m,n)) = (1/dx)*(Uo(Mpos(m,n)) - Uo(Mpos(m,n)-n)); DUy(Mpos(m,n)) = (1/dy)*(Uo(Mpos(m,n)) - Uo(Mpos(m,n)-1));
    DVx(Mpos(m,n)) = (1/dx)*(Vo(Mpos(m,n)) - Vo(Mpos(m,n)-n)); DVy(Mpos(m,n)) = (1/dy)*(Vo(Mpos(m,n)) - Vo(Mpos(m,n)-1));
    
    for i = boundupd
        if Uo(i) >= 0
            DUx(i) = (1/dx)*(Uo(i) - Uo(i-n)); DUy(i) = (1/dy)*(Uo(i+1) - Uo(i));
            DVx(i) = (1/dx)*(Vo(i) - Vo(i-n)); DVy(i) = (1/dy)*(Vo(i+1) - Vo(i));
        else
            DUx(i) = (1/dx)*(Uo(i+n) - Uo(i)); DUy(i) = (1/dy)*(Uo(i+1) - Uo(i));
            DVx(i) = (1/dx)*(Vo(i+n) - Vo(i)); DVy(i) = (1/dy)*(Vo(i+1) - Vo(i));
        end
    end
    for i = bounddownd
        if Uo(i) >= 0
            DUx(i) = (1/dx)*(Uo(i) - Uo(i-n)); DUy(i) = (1/dy)*(Uo(i) - Uo(i-1));
            DVx(i) = (1/dx)*(Vo(i) - Vo(i-n)); DVy(i) = (1/dy)*(Vo(i) - Vo(i-1));
        else
            DUx(i) = (1/dx)*(Uo(i+n) - Uo(i)); DUy(i) = (1/dy)*(Uo(i) - Uo(i-1));
            DVx(i) = (1/dx)*(Vo(i+n) - Vo(i)); DVy(i) = (1/dy)*(Vo(i) - Vo(i-1));
        end
    end     
    for i = boundrightd
        if Vo(i) >= 0
            DUx(i) = (1/dx)*(Uo(i) - Uo(i-n)); DUy(i) = (1/dy)*(Uo(i) - Uo(i-1));
            DVx(i) = (1/dx)*(Vo(i) - Vo(i-n)); DVy(i) = (1/dy)*(Vo(i) - Vo(i-1));
        else
            DUx(i) = (1/dx)*(Uo(i) - Uo(i-n)); DUy(i) = (1/dy)*(Uo(i+1) - Uo(i));
            DVx(i) = (1/dx)*(Vo(i) - Vo(i-n)); DVy(i) = (1/dy)*(Vo(i+1) - Vo(i));
        end
    end    
    for i = boundleftd
        if Vo(i) >= 0
            DUx(i) = (1/dx)*(Uo(i+n) - Uo(i)); DUy(i) = (1/dy)*(Uo(i) - Uo(i-1));
            DVx(i) = (1/dx)*(Vo(i+n) - Vo(i)); DVy(i) = (1/dy)*(Vo(i) - Vo(i-1));
        else
            DUx(i) = (1/dx)*(Uo(i+n) - Uo(i)); DUy(i) = (1/dy)*(Uo(i+1) - Uo(i));
            DVx(i) = (1/dx)*(Vo(i+n) - Vo(i)); DVy(i) = (1/dy)*(Vo(i+1) - Vo(i));
        end
    end
%     for i = bound
%         if Uo(i) >= 0 && Vo(i) >= 0
%             DUx(i) = (1/dx)*(Uo(i) - Uo(i-n)); DUy(i) = (1/dy)*(Uo(i) - Uo(i-1));
%             DVx(i) = (1/dx)*(Vo(i) - Vo(i-n)); DV2(i) = (1/dy)*(Vo(i) - Vo(i-1));
%         elseif Uo(i) > 0 && Vo(i) < 0
%             DUx(i) = (1/dx)*(Uo(i) - Uo(i-n)); DUy(i) = (1/dy)*(Uo(i+1) - Uo(i));
%             DVx(i) = (1/dx)*(Vo(i) - Vo(i-n)); DV2(i) = (1/dy)*(Vo(i+1) - Vo(i));
%         elseif Uo(i) < 0 && Vo(i) > 0
%             DUx(i) = (1/dx)*(Uo(i+n) - Uo(i)); DUy(i) = (1/dy)*(Uo(i) - Uo(i-1));
%             DVx(i) = (1/dx)*(Vo(i+n) - Vo(i)); DV2(i) = (1/dy)*(Vo(i) - Vo(i-1));
%         else
%             DUx(i) = (1/dx)*(Uo(i+n) - Uo(i)); DUy(i) = (1/dy)*(Uo(i+1) - Uo(i));
%             DVx(i) = (1/dx)*(Vo(i+n) - Vo(i)); DV2(i) = (1/dy)*(Vo(i+1) - Vo(i));
%         end
%     end
    for i = boundintd
        if Uo(i) >= 0 && Vo(i) >= 0
            DUx(i) = (1/dx)*(Uo(i) - Uo(i-n)); DUy(i) = (1/dy)*(Uo(i) - Uo(i-1));
            DVx(i) = (1/dx)*(Vo(i) - Vo(i-n)); DVy(i) = (1/dy)*(Vo(i) - Vo(i-1));
        elseif Uo(i) > 0 && Vo(i) < 0
            DUx(i) = (1/dx)*(Uo(i) - Uo(i-n)); DUy(i) = (1/dy)*(Uo(i+1) - Uo(i));
            DVx(i) = (1/dx)*(Vo(i) - Vo(i-n)); DVy(i) = (1/dy)*(Vo(i+1) - Vo(i));
        elseif Uo(i) < 0 && Vo(i) > 0
            DUx(i) = (1/dx)*(Uo(i+n) - Uo(i)); DUy(i) = (1/dy)*(Uo(i) - Uo(i-1));
            DVx(i) = (1/dx)*(Vo(i+n) - Vo(i)); DVy(i) = (1/dy)*(Vo(i) - Vo(i-1));
        else
            DUx(i) = (1/dx)*(Uo(i+n) - Uo(i)); DUy(i) = (1/dy)*(Uo(i+1) - Uo(i));
            DVx(i) = (1/dx)*(Vo(i+n) - Vo(i)); DVy(i) = (1/dy)*(Vo(i+1) - Vo(i));
        end
    end
    for i = 3:m-2
        for j = 3:n-2
            if Uo(i,j) >= 0 && Vo(i,j) > 0
                DUx(i,j) = (1/(2*dx))*(3*Uo(i,j) - 4*Uo(i,j-1) + Uo(i,j-2)); DUy(i,j) = (1/(2*dy))*(3*Uo(i,j) - 4*Uo(i-1,j) + Uo(i-2,j));
                DVx(i,j) = (1/(2*dx))*(3*Vo(i,j) - 4*Vo(i,j-1) + Vo(i,j-2)); DVy(i,j) = (1/(2*dy))*(3*Vo(i,j) - 4*Vo(i-1,j) + Vo(i-2,j));
            elseif Uo(i,j) >= 0 && Vo(i,j) <= 0
                DUx(i,j) = (1/(2*dx))*(3*Uo(i,j) - 4*Uo(i,j-1) + Uo(i,j-2)); DUy(i,j) = (1/(2*dy))*(-Uo(i+2,j) + 4*Uo(i+1,j) - 3*Uo(i,j));
                DVx(i,j) = (1/(2*dx))*(3*Vo(i,j) - 4*Vo(i,j-1) + Vo(i,j-2)); DVy(i,j) = (1/(2*dy))*(-Vo(i+2,j) + 4*Vo(i+1,j) - 3*Vo(i,j));
            elseif Uo(i,j) < 0 && Vo(i,j) > 0
                DUx(i,j) = (1/(2*dx))*(-Uo(i,j+2) + 4*Uo(i,j+1) - 3*Uo(i,j)); DUy(i,j) = (1/(2*dy))*(3*Uo(i,j) - 4*Uo(i-1,j) + Uo(i-2,j));
                DVx(i,j) = (1/(2*dx))*(-Vo(i,j+2) + 4*Vo(i,j+1) - 3*Vo(i,j)); DVy(i,j) = (1/(2*dy))*(3*Vo(i,j) - 4*Vo(i-1,j) + Vo(i-2,j));
            elseif Uo(i,j) < 0 && Vo(i,j) <= 0
                DUx(i,j) = (1/(2*dx))*(-Uo(i,j+2) + 4*Uo(i,j+1) - 3*Uo(i,j)); DUy(i,j) = (1/(2*dy))*(-Uo(i+2,j) + 4*Uo(i+1,j) - 3*Uo(i,j));
                DVx(i,j) = (1/(2*dx))*(-Vo(i,j+2) + 4*Vo(i,j+1) - 3*Vo(i,j)); DVy(i,j) = (1/(2*dy))*(-Vo(i+2,j) + 4*Vo(i+1,j) - 3*Vo(i,j));
            end
        end
    end
end