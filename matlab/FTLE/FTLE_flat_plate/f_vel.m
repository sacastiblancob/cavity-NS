function [U,V] = f_vel(X,Y,t)
% Function for velocity in X direction (U) and velocity in Y direction (V)
%
%   Based on: Exact theory of material spike formation in flow separation.
%               Serra, Vétel, Haller, 2018

Rc=1;       %Cylinder radius
yc=2;       %CYlinder center location
om=20;    %Omega parameter
Uo=0;       %U0 parameter
Uw = -Uo;   %Velocity of the wall
bet=0;      %beta parameter
wc = 1;     %omega sub c, parameter

a = (Rc + yc - sqrt(yc^2-Rc^2))/(Rc + yc + sqrt(yc^2-Rc^2));
mu = 1/(sqrt(yc^2-Rc^2));
ga = (a/(a^2-1))*(-(Uw/(2*log(a))) + (2*om*a^2)/(mu*(a^2-1)^2));

X = X - Uw*t - (bet/wc)*cos(wc*t);
eta = X + 1i*Y;
phi = (1+1i*mu*eta)./(1i+mu*eta);

UV= -(Uw/(2*log(a)))*(2*log(abs(phi)/a) + (mu./(2*phi)).*(conj(eta) - eta).*((phi-1i).^2)) + ...
    (ga*(phi-1i).^2).*((1i*mu*conj(eta)/2).*((a./(phi.^2))+(1/a)) - (1i./phi)*(a + 1/a) + (1/(2*a))*((a^2)./(phi.^2) - 1)) +...
    ga*(a+(1/a)+1i*((a./conj(phi)) - (conj(phi)./a)));

U = real(UV);
V = -imag(UV);

end