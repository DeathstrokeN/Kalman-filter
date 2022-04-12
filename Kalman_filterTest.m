%%
%Kalman filter performances test
%to use the code you need SystemIdentification toolbox :D


%%
clc; clear all;
Ts = 2;
%%
%model
A=0.5;
B=100;
C=1;
D=0;
%%
%Augmanted model
Aa=[A 0; 0 1];
Ba=[B;0];
Ca=[C 1];
N = 50000; % number of points
u=idinput(N,'PRBS',[0 0.5]); % input signal: series of
xp = inv(eye(length(A))-A)*B*u(1); % initial state vector of the plant %u(1)=-1
x_ed = [0;0]; % initial estimates
x_eb = 0;
%Noises
q = 0.5.*randn(1,N) + rand;
p = 0.8.*randn(1,N) + rand;
v = 0.6.*randn(1,N);
w = 0.2.*randn(1,N);
%gain statique du filtre
k = .6;
ka = [k; 0.5];
% Matrix to record the results
% Initialization of matrices to record data
est_errb = zeros(1, N);
est_errd = zeros(1, N);
y_errd = zeros(1, N);
y_errb = zeros(1, N);
%%
% Simulation loop 
for i = 1 : N
    %%
    %d.b
    y = C*xp + v(i) + q(i);
    x_eb = A*x_eb+B*u(i)+k*(y-C*x_eb);
    xp = A*xp + B*u(i) + w(i) + p(i);
    y_errb(i) = y - C*x_eb;
    est_errb(i) = xp - x_eb;
    %%
    %4.d
    x_ed = Aa*x_ed+Ba*u(i)+ka*(y-Ca*x_ed);
    xp = A*xp + B*u(i) + w(i) + p(i);
    y_errd(i) = y - Ca*x_ed;
    est_errd(i) = xp - x_ed(1);
end

display ..
mean(y_errd)
mean(est_errd)


display ..
mean(y_errb)
mean(est_errb)