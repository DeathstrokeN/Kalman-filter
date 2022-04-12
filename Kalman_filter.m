%%
%Simple design of a kalman filter (MIMO sys 2x2)
%with study of :
%               the convergence speed of the estimation errors,
%               the final value of the filter gain, and
%               the final value of the estimated covariance matrix.
% considering the effect of initialization of the covariance matrix

%to use the code you need SystemIdentification toolbox :D

%%
clc; clear all;
Ts = 2;
A=[1.5 -0.6; 1 0];
B=[1; 0];
C=[-0.05 0.5];
D=0;
N = 50000; % number of points
xp = inv(eye(length(A))-A)*B*-1; % initial state vector of the plant %u(1)=-1
x_e = [0; 0]; % initial estimates
Wp = [0.5 0; 0 0.5]; % covariance matrix of the 
 % true (plant) process noise
Vp = 0.8; % covariance matrix of the true (plant) 
 % measurement noise
W = Wp; % estimation of Wp, used in the filter design
V = Vp; % estimation of Vp, used in the filter design
u=idinput(N,'PRBS',[0 0.5]); % input signal: series of
 % steps (pseudo-random
 % binary sequence)
P = 1000*eye(2); % Initialization of P
% Matrix to record the results
% Initialization of matrices to record data
xprocess = zeros(2, N); 
xestimated = zeros(2, N);
est_err = zeros(2, N);
y_err = zeros(1, N);
Pdiag1 = zeros(1, N);
Pdiag2 = zeros(1, N);
%%
%Conception du filtre
% Gkf=eye(2);
% Hkf=zeros(1,2);
% NN=zeros(2,1);
% Gmodel=ss(A,[B Gkf],C,[D Hkf],Ts);
% [KalmanFilter, K, P]=kalman(Gmodel,W,V,NN,'delayed');

%%
% Simulation loop 
for i = 1 : N
    K=A*P*C'*inv(C*P*C'+V);
    v = sqrt(.8)*randn(1,1); % measurement noise
    y = C*xp + v; % plant output
    x_e = A*x_e+B*u(i)+K*(y-C*x_e); % estimation of the states;
    w = sqrt(0.5)*randn(2,1); % process noise
    xp = A*xp + B*u(i) + w; % plant state update
    P=A*(P-P*C'*inv(C*P*C'+V)*C*P)*A'+W;
    
    %collecting data
    xprocess(:,i)=x_e; 
    xestimated(:,i)=xp;
    est_err(:,i) = xp - x_e; 
    y_err(1,i) = y - C*x_e;
    Pdiag1(1,i) = P(1,1);
    Pdiag2(1,i) = P(2,2);
end
%%
%Comparaison*********
KFinal_value = K
%variance de l'erreur d'estimation de x1
var_e_x1=var(est_err(1,:))
p1=P(1,1)
%variance de l'erreur d'estimation de x1
var_e_x2=var(est_err(2,:))
p2=P(2,2)
%*********
mean(y_err)
mean(est_err(1,:))
%%
figure(1);
t=(0:N-1)*Ts;
subplot(211);
plot(t,xprocess(2,:),'k',t,xestimated(2,:),'r');
legend('Process states','Estimated states','south east')
title('comparaison de x2');
subplot(212);
plot(t,xprocess(1,:),'k',t,xestimated(1,:),'r');
title('comparaison de x1');
xlabel('Temps [s]');

figure(2);
plot(t,est_err(2,:),'k',t,est_err(1,:),'r');
legend('on x1','on x2')
title('Estimation error');
xlabel('Temps [s]');


figure(4);
plot(t,y_err);
title('Erreur de prédiction');
xlabel('Temps [s]');


figure(5);
plot(t,Pdiag1,t,Pdiag2);
title('diagonales de P');
xlabel('Temps [s]');