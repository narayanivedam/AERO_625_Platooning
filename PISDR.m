close all;
clear all;
% PI-SDR 
tau1 = 0.1;
tau2 = 0.1;
T = 0.01;
h = 0.01;
t_final = 100;
nframes = t_final/h + 1;
% Continuous time 
A = [0,1,0,0,0;0,-1/tau1,0,0,0;0,0,0,1,0;0,0,0,0,1;0,0,0,0,-1/tau2];
B = [0,0;1/tau1,0;0,0;0,0;0,1/tau2];
C = [1,0,0,0,0;0,0,1,0,0]; % NZSP 
D = 0;
G = [0;1.6;0;0;1.6]; % Owing to wind (Assumed to be a constant)
A_alter = [0,1,0,0,0,0,0;0,-1/tau1,0,0,0,0,0;0,0,0,1,0,0,0;0,0,0,0,1,0,0;0,0,0,0,-1/tau2,0,0];
C_alter = [1,0,0,0,0,0,0;0,0,1,0,0,0,0];
A_New = [A_alter;C_alter];
B_New = [B;0,0;0,0];
sys = ss(A,B,C,D);% State space form
dt_sys = c2d(sys, h, 'zoh');% Continuous to Discrete time system
Phi = dt_sys.a;
Gamma = dt_sys.b;
H = dt_sys.c;
D = dt_sys.d;
GT= G*h;
I = eye(5);
Lower = [H*T,eye(2)];
Upper = [Phi,zeros(5,2)];
Phi_New = [Upper;Lower];
Gamma_New = [Gamma;0,0;0,0];
G_New = [GT;0;0];
% Gains K
Q=[2000 0 0 0 0 0 0 ;0 1 0 0 0 0 0; 0 0 2000 0 0 0 0; 0 0 0 1 0 0 0 ; 0 0 0 0 1 0 0;0 0 0 0 0 150 0; 0 0 0 0 0 0 300 ];
R=[50 0 ; 0 100];
[K,Q_hat,R_hat,M,P,E]=lqrdjv(A_New,B_New,Q,R,T);
% Initial conditions
x0 = [60,0,100,0,0,60,100];
x_ref = [0;0;0;0;0;0;0];
x_ss = [40;0;70;0;0;40;70];% Steady State values of states
x_k = x0'-x_ss;%Change of co-ordinates
x_k1 = x_k;
t_store = [0:h:t_final];
y_store = zeros(2,nframes);
u_store = zeros(2,nframes);
for k = 1:(nframes) 
    if (mod(t_store(k),T) == 0)
        u = K*(x_ref - x_k);
    end
    x_k1 = Phi_New*x_k + Gamma_New*u;
    y = x_k1(6:7);
    y_store(:,k) = y;
    u_store(:,k) = u;
    x_k = x_k1;
end
figure(2);
subplot(2,1,1);
plot(t_store, y_store(1:2,:));
subplot(2,1,2)
plot(t_store, u_store(1:2,:));
% Frequency domain
% H = tf (sys);
% sigma (H) ;