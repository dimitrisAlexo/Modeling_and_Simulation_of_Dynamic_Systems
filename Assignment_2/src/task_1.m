% On-line estimation of unknown parameters using Gradient Descent Method
clear vars;
close all;
clc;
tic;

% Parameters
global a b a_m gamma;
a = 3;
b = 0.5;
a_m = 4;
gamma = 50;

% Gradient Descent Method Estimator
tspan = 0:0.01:30;
y0 = [0 0 0 0 0];
[t,y] = ode45(@(t,y) func(t,y), tspan, y0);    

x = y(:,1);
phi = [y(:,2) y(:,3)];
theta_est = [y(:,4), y(:,5)];

% Estimation of x
x_est = theta_est(:,1).*phi(:,1) + theta_est(:,2).*phi(:,2);

% Estimation of system parameters
a_est = a_m - theta_est(:,1);
b_est = theta_est(:,2);

% Error e
error = x - x_est;

% Lyapunov
V = (1/2)*(a_est-a).^2 + (1/2)*(b_est-b).^2;
Vdot = - gamma*error.^2;

toc;

% Plotting x, x_est, their difference, a_est, and b_est

% Create a figure
fig = figure('Color', 'w');

% Plot x, x_est, and their difference
subplot(2,1,1)
plot(t, x, 'b', 'LineWidth', 1.5);
hold on
plot(t, x_est, 'r', 'LineWidth', 1.5);
plot(t, x - x_est, 'g', 'LineWidth', 1.5);
hold off
xlabel('Time');
ylabel('x, x_{est}, x - x_{est}');
title('Estimation of x');
legend('x', 'x_{est}', 'x - x_{est}');
grid on

% Plot a_est and b_est
subplot(2,1,2)
plot(t, a_est, 'm', 'LineWidth', 1.5);
hold on
plot(t, b_est, 'c', 'LineWidth', 1.5);
hold off
xlabel('Time');
ylabel('a_{est}, b_{est}');
title('Estimation of System Parameters');
legend('a_{est}', 'b_{est}');
grid on

% Set figure properties for better aesthetics
set(fig, 'Position', [100, 100, 800, 600]); % Set figure size
set(fig, 'DefaultLineLineWidth', 1.5); % Set default line width
set(fig, 'DefaultAxesLineWidth', 1.5); % Set default axes line width
set(fig, 'DefaultAxesFontSize', 12); % Set default font size

function dydt = func(t,y)

    global a b a_m gamma;

    dydt = zeros(5,1);
    
    x = y(1);
    phi = [y(2) y(3)];
    theta_est = [y(4) y(5)];
    
    %u = 10;  % Part a
    u = 10*sin(3*t);  % Part b
    
    error = x - (theta_est(1)*phi(1) + theta_est(2)*phi(2));
    
    dydt(1) = - a*x + b*u;
    dydt(2) = - a_m*phi(1) + x;
    dydt(3) = - a_m*phi(2) + u;
    dydt(4) = gamma*error*phi(1);
    dydt(5) = gamma*error*phi(2);
    
end