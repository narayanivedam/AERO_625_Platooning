clear all;
close all;
% Initialization
% For an electric powertrain, the propulsion lag time constant is 0.1s
tau1 = 0.1;
tau2 = 0.1;
% For an ICE, the propulsion lag time constant is 0.5s
% tau1 = 0.5;
% tau2 = 0.5;
h = 0.01;
% Sampling period
T = 0.01;
% Controller Type
c_type = 0;
% Total time
t_final = 10;
% Number of frames
nframes = t_final/h + 1;
% Continuous time 
A = [0,1,0,0,0;0,-1/tau1,0,0,0;0,0,0,1,0;0,0,0,0,1;0,0,0,0,-1/tau2];
B = [0,0;1/tau1,0;0,0;0,0;0,1/tau2];
C = [1,0,0,0,0;0,0,1,0,0]; % NZSP 
D = 0;
sys = ss(A,B,C,D);% State space form
dt_sys = c2d(sys, h, 'zoh');% Continuous to Discrete time system
Phi = dt_sys.a;
Gamma = dt_sys.b;
H = dt_sys.c;
D = dt_sys.d;
I = eye(7);
% Augmenting u to the state vector:
% Initial conditions
x0 = [20,0,70,0,0,0,0];
x_k = x0';
x_k1 = x_k;
t_store = [0:h:t_final];
y_store = zeros(2,nframes);
u_store = zeros(2,nframes);
% New matrices
A_alter = [A,B];
A_New = [A_alter;zeros(2,7)];
B_New = [zeros(5,2);eye(2)];
C_New = [1,0,0,0,0,0,0;0,0,1,0,0,0,0];
sys_New = ss(A_New,B_New,C_New,D);% State space form
dt_sys_New = c2d(sys_New, h, 'zoh');% Continuous to Discrete time system
Phi_New = dt_sys_New.a;
Gamma_New = dt_sys_New.b;
H_New = dt_sys_New.c;
D_New = dt_sys_New.d;
G = [0;1.6;0;0;1.6];
GT= G*h;
G_New = [GT;0;0];
% Controllability check
cc=ctrb(A_New,B_New);% Is controllable
% QPM Inverse
[Pi_12,Pi_22]=QPMCALC(Phi_New-I,Gamma_New,H_New,D_New);
% Gains K
Q=[1000 0 0 0 0 0 0 ;0 10 0 0 0 0 0 ; 0 0 100 0 0 0 0; 0 0 0 1 0 0 0 ; 0 0 0 0 10 0 0 ; 0 0 0 0 0 10 0; 0 0 0 0 0 0 50];
R=[50 0 ; 0 100];
[K,Q_hat,R_hat,M,P,E]=lqrdjv(A_New,B_New,Q,R,T);
% Command Inputs
ym1 = 40;
ym2 = 70;
ym=[ym1;ym2];
% Initial conditions
x0 = [20; 0; 4 ; 0 ;0 ;0;0];
x_k = x0;
x_k1 = x_k;
t_store = [0:h:t_final];
y_store = zeros(2,nframes);
u_store = zeros(2,nframes);
ym_store = zeros(2,nframes);
for k = 1:(nframes) 
    if (mod(t_store(k),T) == 0)
        u = (Pi_22+K*Pi_12)*ym - K*x_k;
    end
    ym=[ym1;ym2]*sin(k*T);
    x_k1 = Phi_New*x_k + Gamma_New*u + G_New;
    y = x_k1(6:7);
    y = H_New*x_k + D_New*u;
    y_store(:,k) = y;
    ym_store(:,k) = ym;
    u_store(:,k) = u;
    x_k = x_k1;
end
figure(2);
subplot(2,1,2)
plot(t_store, u_store(1:2,:));
subplot(2,1,1);
plot(t_store,ym_store(1:2,:),'o');
hold on;
plot(t_store, y_store(1:2,:));

