clc
%% Create the 3-D Simulation Space
view(135,20)
AL = .2;
grid on
%% Redefine the Parameters
L1 = 0.066; %(Meters)
L2 = 0.146; %(Meters)
s = 8;
theta = 0;
phi = 0;
c = [0 0 0];
Xh = [0; L1]';
Yh = [0; 0]';
Zh = [0; 0]';
Xv = [Xh(2); L1]';
Yv = [Yh(2); 0]';
Zv = [Zh(2); -L2]';
hold on
Harm  = fill3(Xh,Yh,Zh,'b'); 
Varm  = fill3(Xv,Yv,Zv,'g'); 
M = scatter3(Xv(2),Yv(2),Zv(2),s);
axis([-AL AL -AL AL -AL AL]);
TXT = title('Time: ');
for t = 1:17:size(simTheta,1)
    TXT2 = sprintf('Time:%.2f',simt(t));
    set(TXT,'String',TXT2);
    phi = simPhi(t);
    theta = -simTheta(t);
    Xva = 0;                
    Yva = L2 * sin(theta);    
    Zva = -L2 * cos(theta);   
    Xvb = Xva * cos(phi) - Yva * sin(phi) + L1 * cos(phi);
    Yvb = Xva * sin(phi) + Yva * cos(phi) + L1 * sin(phi);
    Zvb = Zva;
    Xh(2) = L1 * cos(phi);
    Yh(2) = L1 * sin(phi);
    Xv = [Xh(2); Xvb]';
    Yv = [Yh(2); Yvb]';
    Zv = [0; Zvb]';
    set(Harm,'XData',Xh);
    set(Harm,'YData',Yh);
    set(Harm,'ZData',Zh);
    set(Varm,'XData',Xv);
    set(Varm,'YData',Yv);
    set(Varm,'ZData',Zv);
    rem(t,30) 
    set(M,'XData',Xv(2));
    set(M,'YData',Yv(2));
    set(M,'ZData',Zv(2));
    drawnow;
end