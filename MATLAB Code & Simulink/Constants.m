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
M = 0.044;      % Sin/Cos Values
kb_p = 4.7940e-04 ;
kb_m = 6.75e-4 ;
%% Equations for Lagrange Equation of the system
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
%% Making the Linearization of The System
% Making Matrix A
A = zeros(4,4);
A(1,2) = 1;
A(2,3) = -(sigma * gamma)/(alpha * beta - gamma^2);
A(3,4) = 1;
A(4,3) = (alpha * sigma)/(alpha * beta - gamma^2);
% Making Matrix B
B = zeros(4,2);
B(2,1) = beta/(alpha * beta - gamma^2);
B(2,2) = -gamma/(alpha * beta - gamma^2);
B(4,1) = -gamma/(alpha * beta - gamma^2);
B(4,2) = alpha/(alpha * beta - gamma^2);
% C matrix
C = [0 0 1 0; 0 0 0 1];
%% Linear of The System Matrixs 
% Making New Matrix A
Ap = zeros(4,4);
Ap(1,2) = 1;
Ap(2,1) = 0;      Ap(2,2) = -B(2,1) * (ke^2/Re + kb_m); 
Ap(2,3) = A(2,3); Ap(2,4) = -B(2,2) * kb_p;
Ap(3,4) = 1;
Ap(4,1) = 0;      Ap(4,2) = -B(4,1) * (ke^2/Re + kb_m); 
Ap(4,3) = A(4,3); Ap(4,4) = -B(4,2) * kb_p;
% Making New Matrix B
Bp = zeros(4,1);
Bp(2) = B(2,1) * ke/Re;  
Bp(4) = B(4,1) * ke/Re;   
%% State Feedback Control
K = place(Ap,Bp,[-5 -4 -2+2j -2-2j]);  
Q = [0.1 0 0 0; 0 0.01 0 0; 0 0 100 0; 0 0 0 10];
R = 10;
[K, ~, E] = lqr(Ap,Bp,Q,R);
% Making Controls & Observability
Control = rank(ctrb(Ap,Bp));
Observ = rank(obsv(Ap,C));
%% Pendulum Arm in The Swing Down Position
% Swing Down Matrix A
Ap2 = zeros(4,4);
Ap2(1,2) = 1;
Ap2(2,1) = 0;      Ap2(2,2) = -B(2,1) * (ke^2/Re + kb_m); 
Ap2(2,3) = A(2,3); Ap2(2,4) = B(2,2) * kb_p;
Ap2(3,4) = 1;
Ap2(4,1) = 0;      Ap2(4,2) = B(4,1) * (ke^2/Re + kb_m); 
Ap2(4,3) = -A(4,3); Ap2(4,4) = -B(4,2) * kb_p;
% Swing Down Matrix B
Bp2 = zeros(4,1);
Bp2(2) = B(2,1) * ke/Re;  
Bp2(4) = -B(4,1) * ke/Re;  

K2 = place(Ap2,Bp2,[-5 -4 -2+2j -2-2j]);
R2 = 1;
Q2=[1 0 0 0; 0 10 0 0; 0 0 1000 0; 0 0 0 10];
[K2, ~, E] = lqr(Ap2,Bp2,Q2,R2);