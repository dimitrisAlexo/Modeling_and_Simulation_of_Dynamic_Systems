% On-line estimation of unknown parameters using Lyapunov Method
clear vars;
close all;
clc;
tic;

% Parameters
global a b gamma_par gamma_mix theta_m h0 f
a = 3;
b = 0.5;
gamma_par = [15 1];
gamma_mix = [15 1];
theta_m = 3;
h0 = 0.5;
f = 2;

% Parallel Configuration
tspan = 0:0.001:20;
y0 = [0 0 0 0];
[t_par,y_par] = ode45(@(t,y) lyapunov_par(t,y), tspan, y0);
[t_par_n,y_par_n] = ode45(@(t,y) lyapunov_par_noise(t,y), tspan, y0);

% Without Noise
x = y_par(:,1);
theta_est = [y_par(:,2), y_par(:,3)];
x_est = y_par(:,4);

% With Noise
x_n = y_par_n(:,1);
theta_est_n = [y_par_n(:,2), y_par_n(:,3)];
x_est_n = y_par_n(:,4);

% Unknown parameters estimation
a_est = theta_est(:,1);
b_est = theta_est(:,2);
a_est_n = theta_est_n(:,1);
b_est_n = theta_est_n(:,2);

% Error e
error = x - x_est;
error_n = x_n - x_est_n;

% Lyapunov without Noise
V = (1/2)*(error.^2 + (1/gamma_par(1))*(a_est-a).^2 ...
    + (1/gamma_par(2))*(b_est-b).^2);
Vdot = - a*error.^2;

% Lyapunov with Noise
V_n = (1/2)*(error_n.^2 + (1/gamma_par(1))*(a_est_n-a).^2 ...
    + (1/gamma_par(2))*(b_est_n-b).^2);
Vdot_n = - theta_m*error_n.^2;



% Mixed Configuration
[t_mix,y_mix] = ode45(@(t,y) lyapunov_mix(t,y), tspan, y0);
[t_mix_n,y_mix_n] = ode45(@(t,y) lyapunov_mix_noise(t,y), tspan, y0);

% Without Noise
x_m = y_mix(:,1);
theta_est_mix = [y_mix(:,2), y_mix(:,3)];
x_est_mix = y_mix(:,4);

% With Noise
x_m_n = y_mix_n(:,1);
theta_est_mix_n = [y_mix_n(:,2), y_mix_n(:,3)];
x_est_mix_n = y_mix_n(:,4);

a_est_mix = theta_est_mix(:,1);
b_est_mix = theta_est_mix(:,2);
a_est_mix_n = theta_est_mix_n(:,1);
b_est_mix_n = theta_est_mix_n(:,2);

% Error e
error_mix = x_m - x_est_mix;
error_mix_n = x_m_n - x_est_mix_n;

% Lyapunov without Noise
V_m = (1/2)*(error_mix.^2 + (1/gamma_mix(1))*(a_est_mix-a).^2 ...
    + (1/gamma_mix(2))*(b_est_mix-b).^2);
Vdot_m = - theta_m*error_mix.^2;

% Lyapunov with Noise
V_m_n = (1/2)*(error_mix_n.^2 + (1/gamma_par(1))*(a_est_mix_n-a).^2 ...
    + (1/gamma_par(2))*(b_est_mix_n-b).^2);
Vdot_m_n = - theta_m*error_mix_n.^2;

toc;

% Create Figure 1 and 2: Parallel configuration
figure(1)
subplot(2,1,1)
plot(t_par, x, 'LineWidth', 1.5)
hold on
plot(t_par, x_est, 'LineWidth', 1.5)
plot(t_par, x-x_est, 'LineWidth', 1.5)
hold off
title('Parallel Configuration without Noise')
xlabel('Index')
ylabel('Value')
legend('x', 'x_{est}', 'x - x_{est}')
grid on

subplot(2,1,2)
plot(t_par, a_est, 'LineWidth', 1.5)
hold on
plot(t_par, b_est, 'LineWidth', 1.5)
hold off
title('Parallel Configuration without Noise')
xlabel('Index')
ylabel('Value')
legend('a_{est}', 'b_{est}')
grid on

figure(2)
subplot(2,1,1)
plot(t_par_n, x_n, 'LineWidth', 1.5)
hold on
plot(t_par_n, x_est_n, 'LineWidth', 1.5)
plot(t_par_n, x_n-x_est_n, 'LineWidth', 1.5)
hold off
title('Parallel Configuration with Noise')
xlabel('Index')
ylabel('Value')
legend('x_n', 'x_{est_n}', 'x_n - x_{est_n}')
grid on

subplot(2,1,2)
plot(t_par_n, a_est_n, 'LineWidth', 1.5)
hold on
plot(t_par_n, b_est_n, 'LineWidth', 1.5)
hold off
title('Parallel Configuration with Noise')
xlabel('Index')
ylabel('Value')
legend('a_{est_n}', 'b_{est_n}')
grid on

% Create Figure 3 and 4: Mixed configuration
figure(3)
subplot(2,1,1)
plot(t_mix, x_m, 'LineWidth', 1.5)
hold on
plot(t_mix, x_est_mix, 'LineWidth', 1.5)
plot(t_mix, x_m-x_est_mix, 'LineWidth', 1.5)
hold off
title('Mixed Configuration without Noise')
xlabel('Index')
ylabel('Value')
legend('x', 'x_{est}', 'x - x_{est}')
grid on

subplot(2,1,2)
plot(t_mix, a_est_mix, 'LineWidth', 1.5)
hold on
plot(t_mix, b_est_mix, 'LineWidth', 1.5)
hold off
title('Mixed Configuration without Noise')
xlabel('Index')
ylabel('Value')
legend('a_{est}', 'b_{est}')
grid on

figure(4)
subplot(2,1,1)
plot(t_mix_n, x_m_n, 'LineWidth', 1.5)
hold on
plot(t_mix_n, x_est_mix_n, 'LineWidth', 1.5)
plot(t_mix_n, x_m_n-x_est_mix_n, 'LineWidth', 1.5)
hold off
title('Mixed Configuration with Noise')
xlabel('Index')
ylabel('Value')
legend('x_n', 'x_{est_n}', 'x_n - x_{est_n}')
grid on

subplot(2,1,2)
plot(t_mix_n, a_est_mix_n, 'LineWidth', 1.5)
hold on
plot(t_mix_n, b_est_mix_n, 'LineWidth', 1.5)
hold off
title('Mixed Configuration with Noise')
xlabel('Index')
ylabel('Value')
legend('a_{est_n}', 'b_{est_n}')
grid on

function dydt = lyapunov_par(t,y)

    dydt = zeros(4,1);
    global a b gamma_par
    
    x = y(1);
    theta_est = [y(2) y(3)];
    x_est = y(4);
    
    u = 10*sin(3*t);
    
    error = x - x_est;
    
    dydt(1) = - a*x + b*u;
    dydt(2) = - gamma_par(1)*error*x_est;
    dydt(3) = gamma_par(2)*error*u;
    dydt(4) = -theta_est(1)*x_est + theta_est(2)*u;
end

function dydt = lyapunov_par_noise(t,y)

    dydt = zeros(size(y));
    global a b gamma_par h0 f
    
    % Noise
    h = h0*sin(2*pi*f*t);
    
    x = y(1);
    theta_est = [y(2) y(3)];
    x_est = y(4);
    
    u = 10*sin(3*t);
    
    error = x + h - x_est;
    
    dydt(1) = - a*x + b*u;
    dydt(2) = - gamma_par(1)*error*x_est;
    dydt(3) = gamma_par(2)*error*u;
    dydt(4) = -theta_est(1)*x_est + theta_est(2)*u;
end

function dydt = lyapunov_mix(t,y)

    dydt = zeros(size(y));
    global a b gamma_mix theta_m
    
    x = y(1);
    theta_est = [y(2) y(3)];
    x_est = y(4);
    
    u = 10*sin(3*t);
    
    error = x - x_est;
    
    dydt(1) = - a*x + b*u;
    dydt(2) = - gamma_mix(1)*error*x;
    dydt(3) = gamma_mix(2)*error*u;
    dydt(4) = -theta_est(1)*x + theta_est(2)*u + theta_m*error;
end

function dydt = lyapunov_mix_noise(t,y)

    dydt = zeros(size(y));
    global a b gamma_mix theta_m h0 f
    
    % Noise
    h = h0*sin(2*pi*f*t);
    
    x = y(1);
    theta_est = [y(2) y(3)];
    x_est = y(4);
    
    u = 10*sin(3*t);
    
    error = x + h - x_est;
    
    dydt(1) = - a*x + b*u;
    dydt(2) = - gamma_mix(1)*error*(x+h);
    dydt(3) = gamma_mix(2)*error*u;
    dydt(4) = -theta_est(1)*(x+h) + theta_est(2)*u + theta_m*error;
end