% Least Squares Method

clearvars
clc
close all

% Define system parameters
m = 10;
b = 0.5;
k = 2.5;

% Define input signal
t = 0:0.1:10;
u = 15*sin(3*t) + 8;

% Define initial conditions
y0 = 0;
ydot0 = 0;

% Define differential equation
f = @(t,y) [y(2); (-b*y(2) - k*y(1) + u(round(t*10)+1))/m];

% Solve differential equation using ode45
tic;
[t,y] = ode45(f, t, [y0 ydot0]);
toc;

% Estimating with Least-Squares method
p1 = -1;
p2 = -1;

N = length(t);
phi = zeros(N,3);   
phi(:,1) = lsim(tf([-1 0],[1 -(p1+p2) p1*p2]),y(:,1),t);
phi(:,2) = lsim(tf(-1,[1 -(p1+p2) p1*p2]),y(:,1),t);
phi(:,3) = lsim(tf(1,[1 -(p1+p2) p1*p2]),u,t);

phiT_phi = phi.'*phi;
YT_phi = y(:,1).'*phi;
theta_0 = YT_phi/phiT_phi;

m_est = 1/theta_0(3);
b_est = m_est*(theta_0(1)-(p1+p2));
k_est = m_est*(theta_0(2)+p1*p2);

y_est = theta_0*phi';

% Print parameter estimates
fprintf('m_est = %f\nb_est = %f\nk_est = %f\n', m_est, b_est, k_est);

% Plot input and output signals
plot(t, u, 'r', t, y(:,1), 'b');
title('Input and output signals');
xlabel('Time (s)');
ylabel('Signal');
legend('Input u', 'Output y');

% Plot estimated output and difference from true output
figure
plot(t, y(:,1), 'r', t, y_est, 'b');
title('Comparison between y and y_{est}');
xlabel('Time (s)');
ylabel('Magnitude');
legend('True output y', 'Estimated output y_{est}');

figure
plot(t, y(:,1) - y_est', 'r');
title('Error in estimation');
xlabel('Time (s)');
ylabel('Magnitude');
legend('Difference');
