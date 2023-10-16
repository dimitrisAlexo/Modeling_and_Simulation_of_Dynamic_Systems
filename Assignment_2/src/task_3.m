% On-line estimation of unknown parameters using Lyapunov Method - Second order system
clear vars;
close all;
clc;
tic;

% Parameters
global A b C gamma
a11 = -0.25;
a12 = 3;
a21 = -5;
a22 = 0;
A = [a11 a12; a21 a22];

C = [10 0; 0 10];

b1 = 0.5;
b2 = 1.5;
b = [b1 b2];

gamma = [15 3];

tspan = 0:0.01:20;
y0 = zeros(1,10);
[t,y] = ode45(@(t,y) lyapunov_func(t,y), tspan, y0);

x = [y(:,1) y(:,2)];
x_est = [y(:,3) y(:,4)];

A_est = [y(:,5) y(:,6); y(:,7) y(:,8)];
b_est = [y(:,9) y(:,10)];

error = [x(1) - x_est(1) x(2) - x_est(2)];

P = lyap(A',eye(2));

toc;

% Plot x and x_est for each element
figure;
subplot(2,1,1);
plot(t, x(:,1), 'b', 'LineWidth', 1.5, 'DisplayName', 'x_1');
hold on;
plot(t, x_est(:,1), 'r', 'LineWidth', 1.5, 'DisplayName', 'x1_{est}');
xlabel('Time');
ylabel('x_1');
legend('Location', 'best');
grid on;
title('Estimation of x_1');
subplot(2,1,2);
plot(t, x(:,2), 'b', 'LineWidth', 1.5, 'DisplayName', 'x_2');
hold on;
plot(t, x_est(:,2), 'r', 'LineWidth', 1.5, 'DisplayName', 'x2_{est}');
xlabel('Time');
ylabel('x_2');
legend('Location', 'best');
grid on;
title('Estimation of x_2');
sgtitle('Estimated States x and x_{est}', 'FontSize', 16);


% Plot the difference between x and x_est for each element
figure;
subplot(2,1,1);
plot(t, x(:,1) - x_est(:,1), 'm', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Error in x_1');
grid on;
title('Error in x_1 Estimation');
subplot(2,1,2);
plot(t, x(:,2) - x_est(:,2), 'm', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Error in x_2');
grid on;
title('Error in x_2 Estimation');
sgtitle('Error in Estimated States x and x_{est}', 'FontSize', 16);


% Plot A for each element
figure;
subplot(2,2,1);
plot(t, y(:,5), 'b', 'LineWidth', 1.5);
xlabel('Time');
ylabel('a_1_1_{est}');
grid on;
title('Estimation of a_1_1');
subplot(2,2,2);
plot(t, y(:,6), 'b', 'LineWidth', 1.5);
xlabel('Time');
ylabel('a_1_2_{est}');
grid on;
title('Estimation of a_1_2');
subplot(2,2,3);
plot(t, y(:,7), 'b', 'LineWidth', 1.5);
xlabel('Time');
ylabel('a_2_1_{est}');
grid on;
title('Estimation of a_2_1');
subplot(2,2,4);
plot(t, y(:,8), 'b', 'LineWidth', 1.5);
xlabel('Time');
ylabel('a_2_2_{est}');
grid on;
title('Estimation of a_2_2');
sgtitle('Estimated Parameters A', 'FontSize', 16);


% Plot b for each element
figure;
subplot(2,1,1);
plot(t, y(:,9), 'r', 'LineWidth', 1.5);
xlabel('Time');
ylabel('b_1_{est}');
grid on;
title('Estimation of b_1');
subplot(2,1,2);
plot(t, y(:,10), 'r', 'LineWidth', 1.5);
xlabel('Time');
ylabel('b_2_{est}');
grid on;
title('Estimation of b_2');
sgtitle('Estimated Parameters b', 'FontSize', 16);


function dydt = lyapunov_func(t,y)

    global A b C gamma
    
    u = 3.5*sin(7.2*t) + 2*sin(11.7*t);
    
    x = [y(1) y(2)];
    x_est = [y(3) y(4)];
    
    a11_est = y(5);
    a12_est = y(6);
    a21_est = y(7);
    a22_est = y(8);
    b1_est = y(9);
    b2_est = y(10);
  
    dydt = zeros(size(y));

    error1 = x(1) - x_est(1);
    error2 = x(2) - x_est(2);

    dydt(1) = A(1,1)*x(1) + A(1,2)*x(2) + b(1)*u;
    dydt(2) = A(2,1)*x(1) + A(2,2)*x(2) + b(2)*u;
    
    dydt(3) = a11_est*x(1) + a12_est*x(2) + b1_est*u + C(1,1)*error1 + C(1,2)*error2;
    dydt(4) = a21_est*x(1) + a22_est*x(2) + b2_est*u + C(2,1)*error1 + C(2,2)*error2;
    
    dydt(5) = gamma(1)*x(1)*error1;
    dydt(6) = gamma(1)*x(2)*error1;
    dydt(7) = gamma(1)*x(1)*error2;
    dydt(8) = gamma(1)*x(2)*error2;
    
    dydt(9) = gamma(2)*u*error1;
    dydt(10) = gamma(2)*u*error2;
end