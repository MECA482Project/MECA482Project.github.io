**By John Grenard & Travis Gubbins**

**Spring Semester 2022**

**MECA 482-01** 

![Furuta Pendulum Picture](https://github.com/MECA482Project/MECA482Project.github.io/blob/main/Pictures%20%26%20GIFs/Picture%20of%20Furuta%20Pendulum.PNG)\

**Figure 1.** Furuta Pendulum

You can use the [editor on GitHub](https://github.com/MECA482Project/MECA482Project.github.io/edit/main/README.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

# The 5 Minute Presentation 

# MECA 482 Project Report 

# MATLAB Code For Creating The Constants

Shown below is the MATLAB code used in creating the mathematical model of the system.

```markdown
clear
clc
%% Define the Parameters of The System
g = 9.81;       % (Meter/Second^2)
m1 = 0.380;     % (Kilograms)
m2 = 0.054;     % (Kilograms)
L1 = 0.066;     % (Meters)
L2 = 0.146;     % (Meters)
J = 3.5256e-4;  % (Kilogram-Meter^2)
Re = 14.5;      % Potential energy
ke = 0.5;       % Kinetic Energy
M = 0.044; 
kb_p = 4.7940e-04 ;
kb_m = 6.75e-4 ;

alpha = J + (M + m1/3 + m2) * L1^2;
beta = (M + m2/3) * L2^2;
gamma = (M + m2/2) * L2*L1;
sigma = (M + m2/2) * g * L2;
%% Define the Parameters of The Simulation
Zn = 3;
Ts = .001;
StepX = 10;
dtDisc = 0.01;
initial_state = pi;
distrub = 12;
disturb = distrub*pi/180;
Reference = [0 0 0 0];
%% Linearization
% A matrix
A = zeros(4,4);
A(1,2) = 1;
A(2,3) = -(sigma * gamma)/(alpha * beta - gamma^2);
A(3,4) = 1;
A(4,3) = (alpha * sigma)/(alpha * beta - gamma^2);
% B matrix
B = zeros(4,2);
B(2,1) = beta/(alpha * beta - gamma^2);
B(2,2) = -gamma/(alpha * beta - gamma^2);
B(4,1) = -gamma/(alpha * beta - gamma^2);
B(4,2) = alpha/(alpha * beta - gamma^2);
% C matrix
C = [0 0 1 0; 0 0 0 1];
%% Pseudo Linear System
% Pseudo A matrix
Ap = zeros(4,4);
Ap(1,2) = 1;
Ap(2,1) = 0;      Ap(2,2) = -B(2,1) * (ke^2/Re + kb_m); 
Ap(2,3) = A(2,3); Ap(2,4) = -B(2,2) * kb_p;
Ap(3,4) = 1;
Ap(4,1) = 0;      Ap(4,2) = -B(4,1) * (ke^2/Re + kb_m); 
Ap(4,3) = A(4,3); Ap(4,4) = -B(4,2) * kb_p;
% Pseudo B matrix
Bp = zeros(4,1);
Bp(2) = B(2,1) * ke/Re;  
Bp(4) = B(4,1) * ke/Re;   
% Controlability and Observability
Control = rank(ctrb(Ap,Bp));
Observ = rank(obsv(Ap,C));
% State Feedback Control
K = place(Ap,Bp,[-5 -4 -2+2j -2-2j]);  
Q = [0.1 0 0 0; 0 0.01 0 0; 0 0 100 0; 0 0 0 10];
R = 10;
[K, ~, E] = lqr(Ap,Bp,Q,R);
%% Down position (pi)
% Pseudo A matrix around pi
Ap2 = zeros(4,4);
Ap2(1,2) = 1;
Ap2(2,1) = 0;      Ap2(2,2) = -B(2,1) * (ke^2/Re + kb_m); 
Ap2(2,3) = A(2,3); Ap2(2,4) = B(2,2) * kb_p;
Ap2(3,4) = 1;
Ap2(4,1) = 0;      Ap2(4,2) = B(4,1) * (ke^2/Re + kb_m); 
Ap2(4,3) = -A(4,3); Ap2(4,4) = -B(4,2) * kb_p;
% Pseudo B matrix around pi
Bp2 = zeros(4,1);
Bp2(2) = B(2,1) * ke/Re;  
Bp2(4) = -B(4,1) * ke/Re;  

K2 = place(Ap2,Bp2,[-5 -4 -2+2j -2-2j]);
R2 = 1;
Q2=[1 0 0 0; 0 10 0 0; 0 0 1000 0; 0 0 0 10];
[K2, ~, E] = lqr(Ap2,Bp2,Q2,R2);
```

# MATLAB Code For Creating The Simulation

Shown below is the MATLAB code used in creating the simulation of the system.

```markdown
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

M = scatter3(Xv(2),Yv(2),Zv(2),s);     %,'filled','MarkerFaceColor','r','MarkerEdgeColor','k');

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
```

# Simulink Model

![Simulink Model](https://github.com/MECA482Project/MECA482Project.github.io/blob/main/Screenshot%20of%20the%20Simulink.PNG)

**Figure 2.** The Simulink Model

# MATLAB Simulation Based Off Our Simulink Model

![MATLAB Simulation Based Off Our Simulink Model](https://github.com/MECA482Project/MECA482Project.github.io/blob/main/482%20Project%20Pendulum%20GIF.gif)

**Figure 3.** The MATLAB Simulation


# References

[1] Norman S. Nise - Control Systems Engineering-Wiley (2015) 7th Edition <br>

[2] Lecture and reference videos from H. Sinan Bank’s Blackboard page <br>

[3] Wikipedia, Furuta Pendulum, Retrieved by Apr., 19, 2022 from https://en.wikipedia.org/wiki/Furuta_pendulum <br>

[4] IEEEXPLORE, Modeling, Simulation, and Construction of a Furuta Pendulum Test-Bed, Retrieved by Apr., 25, 2022 from <br>
https://ieeexplore.ieee.org/document/7086928 <br>

[5] Control System Tutorials for MATLAB and Simulink, Retrieved by Apr, 22, 2022 from <br>
http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling <br>

[6] GitHub Docs, Creating and highlighting code blocks, Retrieved by May, 02, 2022 from <br>
https://docs.github.com/en/get-started/writing-on-github <br>

[7] Quanser, QUBE – Servo 2, Retrieved by May, 02, 2022 from <br>
https://www.quanser.com/products/qube-servo-2/ <br>

[8] Screencast-O-Matic, Sites used on May, 19, 2022 from <br> 
https://screencast-o-matic.com/ <br>
