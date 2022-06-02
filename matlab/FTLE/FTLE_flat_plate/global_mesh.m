% Nxg = 10;
% Nyg = 3;
% 
% LXG = ((xmax-xmin)/Nxg)*ones(1,Nxg);
% LYG = ((ymax-ymin)/Nyg)*ones(Nyg,1);
% 
% Nx = 4;
% Ny = 4;
% ep = (0.5 + 0.5*JacobiGL(0,0,Nx))';
% et = (0.5 + 0.5*JacobiGL(0,0,Ny));
% 
% %Coordinates of X-direction element boundaries
% XBH = zeros(1,Nx*Nxg+1);
% YBH = zeros(Nyg+1,1);
% 
% %Coordinates of Y-direction element boundaries
% XBV = zeros(1,Nxg+1);
% YBV = zeros(Ny*Nyg+1,1);
% 
% for i=1:Nxg
%     XBH((i-1)*Nx+1:i*Nx+1) = LXG(i)*ep + (xmin + sum(LXG(1:i-1)));
%     XBV(i) = (xmin + sum(LXG(1:i-1)));
% end
% for i=1:Nyg
%     YBV((i-1)*Ny+1:i*Ny+1) = LYG(i)*et + (ymin + sum(LYG(1:i-1)));
%     YBH(i) = (ymin + sum(LYG(1:i-1)));
% end
% XBV(Nxg+1) = xmin + sum(LXG);
% YBH(Nyg+1) = ymin + sum(LYG);
% 
% 
% XBH = kron(ones(Nyg+1,1),XBH);
% YBH = kron(ones(1,Nx*Nxg+1),YBH);
% 
% YBV = kron(ones(1,Nxg+1),YBV);
% XBV = kron(ones(Ny*Nyg+1,1),XBV);
% 
% plot(XBH,YBH,'k')
% hold on
% plot(XBV',YBV','b')
% axis equal