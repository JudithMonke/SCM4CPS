clc
clear all

Fc = 60; %frequency 60Hz
Av =1;  % amplitude
Ts = 0.0005; % 1/0.0005 = 2 kHz  campionamento
StopTime = 0.1;             % seconds
t = (0:Ts:StopTime-Ts)';     % seconds
omega = 2*pi*Fc; % angular frequency constant over time
phi = 0;

%% 
% $$V\left(t\right)=A_{v\;} \cos \left(\omega t+\phi \right)=A_V \cos \left(\omega 
% t\right)\;\cos \left(\phi \right)-A_V \;\sin \left(\omega t\right)\;\sin \left(\phi 
% \right)$$
% 
% $\left\lbrack \begin{array}{c}x_1 \left(t+1\right)\\x_2 \left(t+1\right)\end{array}\right\rbrack 
% =\left\lbrack \begin{array}{cc}1 & 0\\0 & 1\end{array}\right\rbrack \left\lbrack 
% \begin{array}{c}x_1 \left(t\right)\\x_2 \left(t\right)\end{array}\right\rbrack 
% +w\left(t\right)$  --> $x\left(t+1\right)=\textrm{Ax}\left(t\right)+w\left(t\right)$     
% dim = n
% 
% $y\left(t\right)=\left\lbrack \begin{array}{cc}\cos \left(\omega t\right) 
% & -\sin \left(\omega t\right)\end{array}\right\rbrack \left\lbrack \begin{array}{c}x_1 
% \left(t\right)\\x_2 \left(t\right)\end{array}\right\rbrack +v\left(t\right)$   
% -->  $y\left(t\right)=C\;x\left(t\right)+v\left(t\right)$     dim = m

A =  [1 0; 0 1]; % valori paper: A = eye(2);
B = [0;0]; % valori paper: B = 0;
C =  [cos(omega*t), -sin(omega*t)];  % frequenza di 60Hz
D = 0;
W = eye(2);
U = 1;
n = 2;
m = 1;
Q = eye(2); %state noise covariance
R = 1;  % measurement noise covariance

%to semplify the plot use 0 noise
y_noise = 0.02 * randn(size(t));
x_noise = 0.0 * randn(size(t));

x1 = 1*cos(phi);     % assume initial phase = 0
x2 = 1*sin(phi);     
x(:,1) = [x1;x2];
for i=2:length(t)
    x(:,i) = A*x(:,i-1)  + [x_noise(i); x_noise(i)] ;    % assume 0 noise on x for semplicity
end

for i=1:length(t)
    y(i) = C(i,:)*x(:,i) + y_noise(i);
end

% Signal without attacks

[y_hat, residue, p0, K,  x_hat_Return, K_steps, P_steps]= kalmanfilter(y, Q, R, t);
residue = y - y_hat;

g0 = compute_g(C, t, P_steps, residue);

euc0 = euclidean_distance(y, y_hat, t);

figure(1);
subplot(311)
plot(t,y, t, y_hat, t, zeros(length(t)), 'yellow');
legA = legend('real','estimated','attack signal');
ylim([-1.5,1.5]);
title('signal without attack');grid on;
figure(1)
subplot(312)
plot(t, g0);
ylim([-0.5,4]);
yline(3,'-.k','Threshold');
title('Chi-square detector');
figure(1)
subplot(313)
plot(t, euc0);
ylim([-0.5,4]);
yline(0.6,'-.k','Threshold');
title('Euclidean Detector');

%% Signal w/ FDIA
% *Theorem. 2* [Sinopoli]: given that $C\upsilon \in \textrm{span}\left(\Gamma 
% \right)$ , exists$y^*$ s.t.  $\Gamma y^* =\textrm{Cv}$. 
% 
% v is the unstable eigenvector for A.
% 
% $\Gamma =\textrm{diag}\left(\gamma_1 ,\ldotp \ldotp \ldotp ,\gamma_m \;\right)$ 
% where m = number of sensors = 1.  $\gamma_i =1\;\textrm{if}\in S_{\textrm{bad}}$. 

lambda = 1; % autovalore di A
v = [0;1]; % autovettore di A
Gamma = 1;

% Construction of ya_0 and ya_1
% C = C(1:2, :);   %bring matric C at instant i=1, i=2 
%% 
% from eq (28) [Sinopoli] we found {y0,...yn-1} , with n=2, as done in the code. 
% 
% $$\Delta e(k+1) = (A-KC(k+1) A) \Delta e(k) - K \Gamma y_a(k+1)   \ \  \ (28)$$
% 
% $$\Delta z(k+1) = C(k+1) A \Delta e(k) + \Gamma y_a(k+1)    \    \    \    
% (29)$$

%from 
% y_prime = y
% y = C x_prime + Gamma ya
% Co_1 * ya = v


% [P,~,~] = idare(A',C',Q,R);
% K = P*C'*inv(C*P*C'+R);
% Gamma = 1;

Co = -[(A - K_steps(:,2) * C(2,:) * A)*K_steps(:,1)*Gamma, K_steps(:,2) * Gamma];

% Co = -[(A - K*C* A)*K*Gamma, K* Gamma];

ya = Co\v;


e_1 = -K_steps(:,1) * Gamma * ya(1);
e_2 = (A - K_steps(:,2) * C(2,:) * A)*e_1 - K_steps(:,2) * Gamma * ya(2) %=v

z_1 =  Gamma * ya(1);
z_2 = C(2,:) * A * e_1 + Gamma * ya(2);


%% 
% $$M = max_{k=0,...n-1} || \Delta z(k) ||    \    \    \                                    
% (33)$$
% 
% _attack sequence:_ $y_a(n+i) = y_a(i) - ( \lambda ^{i+1} / M) y^{*} \ \ \     
% i=0,1,...       \ \ \ (34)$

M = max(norm(z_1),norm(z_2));
ya = ya/M;

e_1 = -K_steps(:,1) * Gamma * ya(1);
e_2 = (A - K_steps(:,2) * C(2,:) * A)*e_1 - K_steps(:,2) * Gamma * ya(2);

e(:,1) = e_1;
e(:,2) = e_2;

z(1) = z_1;
z(2) = z_2;

% z(i) = C(i,:) * A * e(:,i-1) + Gamma * ya(2);

for i=3:length(t)
    e(:,i) = e(:,i-1) + (lambda^(i+1) / M * v);
end

%% 
% $$y'(t) = C(t)x'(t) + v(t) + \Gamma y_a(t)$$


yA = y;
for i=3:length(t)
    yA(i) = C(i,:)*x(:,i) + t(i)*10;
end

[yA_hat, residue_A, ~, ~,  x_hat_Return_A, ~, P_steps_A]= kalmanfilter(yA, Q, R,t);
residue_A = yA - yA_hat;

for i=1:length(t)
yA_hat(i) = C(i,:)*x_hat_Return_A(:,i);
end
gA = compute_g(C, t, P_steps_A, residue_A);
euc_A = euclidean_distance(y, yA_hat, t);

delta_z = residue_A - residue;

delta_x1 = x_hat_Return_A(1,:) - x_hat_Return(1,:);
delta_x2 = x_hat_Return_A(2,:) - x_hat_Return(2,:);
for i=1:length(t)
    norm_x1(i) = norm(delta_x1(i));
    norm_x2(i) = norm(delta_x2(i));
    norm_z(i) = norm(delta_z(i));
end


figure(3)
subplot(511)
plot(t, y, 'blue', t, yA_hat, 'red', t, t*10, 'yellow');
ylim([-2,2]);
title('FDIA');grid on;

subplot(512)
plot(t,gA, 'LineWidth',2);
yline(0,'-');
yline(3,'-.k','Threshold');
ylim([0, 3.5])
title('Chi-square detector');
grid on;

subplot(513)
plot(t, euc_A);
yline(0,'-');
yline(0.6,'-.k','Threshold');
title('Euclidean detector');
grid on;

% % Check theorem 
% figure()
% subplot(211)
% hold on;
% plot(norm_x1); %converge
% plot(norm_x2);
% ylim([-1,2]);
% grid on;
% title('|| \Delta x_1 || and || \Delta x_2 ||');
% hold off;

subplot(514)
plot(t,e);
title('|| \Delta e ||');
grid on;

subplot(515)
plot(norm_z);
yline(1,'-.k', 'Limit for \Delta z');
ylim([-1,2]);
title('|| \Delta z ||');
grid on;

norm_z(i) <= 1

%% ATTACK ON STATE 
% in this way x' goes away from the nominal x while remaining undetectable.

% y = C*x;
% ya = y + a;
% a = C*c;
   
x_new1 = x(1,:);
x_new2 = x(2,:);

x_new(1,:) = 1*x_new1;
x_new(2,:) = 1*x_new2;
for i=2:length(t)
    x_new(1,i) = 1*x_new(1,i-1) + 0.005;
    x_new(2,i) = 1*x_new(2,i-1) ;
end


% C*x_prime + ya = y --> ya = y - C*x_prime
% 
% for i=1:length(t)
%     y_a(i) = y(i) - C(i,:)*x_new(:,i);
% end
% 
% for i=1:length(t)
%     y_new(i)= C(i,:)*x_new(:,i) + y_a(i);
% end

for i=1:length(t)
    y_new(i)= C(i,:)*x_new(:,i);
end

[y_hat_new, ~ , p0_a, K_a,  x_hat_Return_new, K_steps_new, P_steps_new]= kalmanfilter(y_new, Q, R, t);

residue_new = y_new - y_hat_new;

delta_r_new = residue_new - residue;

euc_new = euclidean_distance(y,y_hat_new , t);
gN = compute_g(C, t, P_steps_new, residue_new);


figure(4)
subplot(511)
plot(t, y, 'blue', t, y_hat_new, 'red', t, x_new, 'yellow');
ylim([-2,2]);
title('FDIA on state attack');grid on;

subplot(512)
plot(t, gN, 'LineWidth',2);
yline(0,'-');
yline(3,'-.k','Threshold');
title('Chi-square detector');
grid on;

subplot(513)
plot(t, euc_new);
yline(0,'-');
yline(0.6,'-.k','Threshold');
title('Euclidean detector');
grid on;

delta_z_new = residue_new - residue;
delta_x1_new = x_new(1,:) - x(1,:);
delta_x2_new = x_new(2,:) - x(2,:);
for i=1:length(t)
    norm_x1_new(i) = norm(delta_x1_new(i));
    norm_x2_new(i) = norm(delta_x2_new(i));
    norm_z_new(i) = norm(delta_z_new(i));
end

% figure(6)
% Check theorem 
subplot(514)
hold on;
plot(norm_x1_new); %converge
plot(norm_x2_new);
ylim([-1,2]);
grid on;
title('|| \Delta x_1 || and || \Delta x_2 ||');
hold off;

subplot(515)
plot(norm_z_new);
yline(1,'-.k', 'Limit for \Delta z');
ylim([-1,2]);
title('|| \Delta z ||');
grid on;

norm_z_new(i) <= 1
%% 
% 
% 
% CASE  $y = H x + a$
% 
% where $a = H c$

a = C * v;

yAA = y;
for i=3:length(t)
    yAA(i) = y(i) + a(i);
end


[yAA_hat, residue_AA, ~, ~,  x_hat_Return_AA, ~, P_steps_AA]= kalmanfilter(yAA, Q, R,t);
residue_AA = yAA - yAA_hat;

for i=1:length(t)
    yAA_hat(i) = C(i,:)*x_hat_Return_AA(:,i);
end
gAA = compute_g(C, t, P_steps_AA, residue_AA);
euc_AA = euclidean_distance(y, yAA_hat, t);

delta_zA = residue_AA - residue;


figure(7)
subplot(411)
plot(t, y, 'blue', t, yAA_hat, 'red', t, a, 'yellow');
ylim([-2,2]);
title('FDIA y = Cx + a');grid on;

subplot(412)
plot(t,gAA, 'LineWidth',2);
yline(0,'-');
yline(3,'-.k','Threshold');
ylim([-1, 3.5])
title('Chi-square detector');
grid on;

subplot(413)
plot(t, euc_AA);
yline(0,'-');
yline(0.6,'-.k','Threshold');
title('Euclidean detector');
grid on;

% for i=1:length(t)
%     error_1(i) = (x(1,i) - x_hat_Return_AA(1,i)) - (x(1,i) - x_hat_Return(1,i));
%     error_2(i) = (x(2,i) - x_hat_Return_AA(2,i)) - (x(2,i) - x_hat_Return(2,i));
% end
% 

delta_x1A = x_hat_Return_AA(1,:) - x_hat_Return(1,:);
delta_x2A = x_hat_Return_AA(2,:) - x_hat_Return(2,:);
for i=1:length(t)
    norm_x1A(i) = norm(delta_x1(i));
    norm_x2A(i) = norm(delta_x2(i));
    norm_zA(i) = norm(delta_zA(i));
end

%residuo 
subplot(414)
plot(t, norm_zA);
yline(1,'-.k', 'Limit for \Delta z');
ylim([-1,2]);
title('|| \Delta z ||');
grid on;

norm_zA(i) <= 1
%% 
% 
% Kalman

    %this Kalman filter has been modified to follow an input siusoidal
    %"signal" with approximatley 60 Hz.
    %Q_in and R_in are tuning parameters.-
    %filtered_sig is the Kalman-filtered output of the signal.
%%
    function [y_hat, residue, p0, K,  x_hat_Return, K_steps, P_steps]= kalmanfilter(y, Q, R,simulation_time)

    A=[1,0;0,1];
    w=2*pi*60;
    C = [cos(w*simulation_time) -sin(w*simulation_time)];

    %initilize values
    x0 = [0;0];   %assume initial state = 0;
    p0 = eye(size(A));
    Q=Q;
    R=R;
    K = 0;

    for i=1:length(simulation_time)
        %1 prediction error covariance P(k|k-1) = A P(k-1) A'+ Q
        pre_error_cov = A*p0*A'+Q;
        
        %2 kalman gain equation k(k) = P(k|k-1) C' (C P(k|k-1) C' + R)^-1
        K = pre_error_cov*C(i,:)'*inv(C(i,:)*pre_error_cov*C(i,:)'+ R);
        
        %3 filter equation x_hat(k) = x_hat(k|k-1) + k(k)( y(k) - C x_hat(k|k-1) )
        x_hat = A*x0 + K*(y(i)-C(i,:)*A*x0);
        
        %4 error covariance P(k) = (I - k(k) C) P(k|k-1)
        error_cov = pre_error_cov-K*C(i,:)*pre_error_cov;
        
        %recreating the sinusoid   y_hat = C x_hat
        y_hat(i)=C(i,:)*x_hat;

        % residuo: r(t) = y(t)-y_hat(t|t-1)
        residue(i) = y(i)-C(i,:)*A*x0;
        x_hat_Return(:,i) = x_hat;
        x0=x_hat;
        p0=error_cov;   
        P_steps(:,:,i)=error_cov;   
        K_steps(:,i) = K;
        
    end
end
%% 
% EUCLIDIAN

function [distance_vector] = euclidean_distance(amplitudes_signal, amplitudes_estimate, time)
    if (size(amplitudes_signal) ~= size(amplitudes_estimate))
        error('size of true signal and size of estimated signal are not equal');
    end
    distance_vector = zeros(size(amplitudes_signal));
    for t=1:length(distance_vector)
        p = [time(t) amplitudes_signal(t)];
        q = [time(t) amplitudes_estimate(t)];
        %distance_vector(t) = sqrt((amplitudes_signal(t) - amplitudes_estimate(t)^2));
        distance_vector(t) = sqrt((p(1) - q(1))^2 + (p(2) - q(2))^2); % formula della distanza tra due punti
        %distance_vector(t) = dtw(p, q, 'euclidean');
    end
end
%% 
% CHI

function [g] = compute_g(C, time, P_steps, residue)
    for i=1:length(time)
        g_prova(i) = C(i,:)*P_steps(:,:,1)*C(i,:)';
    end
    for i=1:length(time)
        g(i) = residue(i) * g_prova(i) * residue(i);
    end
end

%% 
% 
% 
% 
% 
% 
% 
%