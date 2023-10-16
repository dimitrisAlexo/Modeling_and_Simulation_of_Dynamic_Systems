% Transfer Table Estimation with Least-Squares Method

clearvars
clc
close all

% Define input/output signal
t = 0:0.000001:5;
N = length(t);
u_1 = 2*sin(4*t);
u_1d = 8*cos(4*t);
u_2 = 4*ones(1, N);
u_2d = zeros(1, N);
[V_R, V_C] = v(t);

% Estimating with Least-Squares method
p1 = 400;
p2 = 400;

phi = zeros(N,6);
phi(:,1) = lsim(tf([-1 0],[1 (p1+p2) p1*p2]),V_C,t);
phi(:,2) = lsim(tf(-1,[1 (p1+p2) p1*p2]),V_C,t);         
phi(:,3) = lsim(tf([1 0],[1 (p1+p2) p1*p2]),u_1,t);
phi(:,4) = lsim(tf(1,[1 (p1+p2) p1*p2]),u_1,t);
phi(:,5) = lsim(tf([1 0],[1 (p1+p2) p1*p2]),u_2,t);
phi(:,6) = lsim(tf(1,[1 (p1+p2) p1*p2]),u_2,t);

phiT_phi = phi.'*phi;
VCT_phi = V_C'.'*phi;
theta_0 = VCT_phi/phiT_phi;
theta_h = theta_0 + [p1+p2 p1*p2 0 0 0 0];

RC_est = 1/((1/3)*(theta_h(1)+theta_h(3)+theta_h(5)));
LC_est = 1/((1/2)*(theta_h(2)+theta_h(6)));

V_C_est = theta_0*phi';

V_R_est = u_1+u_2-V_C_est;

% Part b - Noise / Ruido

VC_b = V_C_est;
VR_b = V_R_est;

for i = 1:1:3
    j = randi(N);
    VC_b(j) = VC_b(j) + randi([1000,10000]);
    j = randi(N);
    VR_b(j) = VR_b(j) + randi([1000,10000]);
end

phi_b = zeros(N,6);
phi_b(:,1) = lsim(tf([-1 0],[1 (p1+p2) p1*p2]),VC_b,t);
phi_b(:,2) = lsim(tf(-1,[1 (p1+p2) p1*p2]),VC_b,t);
phi_b(:,3) = lsim(tf([1 0],[1 (p1+p2) p1*p2]),u_1,t);
phi_b(:,4) = lsim(tf(1,[1 (p1+p2) p1*p2]),u_1,t);
phi_b(:,5) = lsim(tf([1 0],[1 (p1+p2) p1*p2]),u_2,t);
phi_b(:,6) = lsim(tf(1,[1 (p1+p2) p1*p2]),u_2,t);

phiT_phi_b = phi_b.'*phi_b;
VCT_phi_b = VC_b'.'*phi_b;
theta_0_b = VCT_phi_b/phiT_phi_b;
theta_h_b = theta_0_b + [p1+p2 p1*p2 0 0 0 0];

RC_est_b = 1/((1/3)*(theta_h_b(1)+theta_h_b(3)+theta_h_b(5)));
LC_est_b = 1/((1/2)*(theta_h_b(2)+theta_h_b(6)));

V_C_est_b = theta_0_b*phi_b';

V_R_est_b = u_1+u_2-V_C_est_b;

% Print parameter estimates
fprintf('RC_est = %f\nLC_est = %f\n', RC_est, LC_est);

% Print circuit transfer table
fprintf("Transfer Function Matrix:\n")
G1 = tf([1 0 1/LC_est],[1 1/RC_est 1/LC_est]);
G2 = tf([1 0 0],[1 1/RC_est 1/LC_est]);
G3 = tf([1/RC_est 0] ,[1 1/RC_est 1/LC_est]);
G4 = tf([1/RC_est 1/LC_est],[1 1/RC_est 1/LC_est]);
G = [G1 G2; G3 G4]

% Plot input and output signals
figure
plot(t, u_1, 'r', t, u_2, 'b');
title('Input signals');
xlabel('Time (s)');
ylabel('Signal');
legend('Input u_{1}', 'Input u_{2}');

figure
plot(t, V_R, 'r', t, V_C, 'b');
title('Output signals');
xlabel('Time (s)');
ylabel('Signal');
legend('Output V_{R}', 'Output V_{C}');

% Plot estimated V_C and difference from true V_C
figure
plot(t, V_C, 'r', t, V_C_est', 'b');
title('Comparison between V_{C} and V_{Cest}');
xlabel('Time (s)');
ylabel('Magnitude');
legend('True output V_{C}', 'Estimated output V_{Cest}');

figure
plot(t, V_C - V_C_est, 'r');
title('Error in estimation of V_{C}');
xlabel('Time (s)');
ylabel('Magnitude');
legend('Difference');

% Plot estimated V_R and difference from true V_R
figure
plot(t, V_R, 'r', t, V_R_est, 'b');
title('Comparison between V_{R} and V_{Rest}');
xlabel('Time (s)');
ylabel('Magnitude');
legend('True output V_{R}', 'Estimated output V_{Rest}');

figure
plot(t, V_R - V_R_est, 'r');
title('Error in estimation of V_{R}');
xlabel('Time (s)');
ylabel('Magnitude');
legend('Difference');

% Part b - plots

% Print parameter estimates
fprintf('RCb_est = %f\nLCb_est = %f\n', RC_est_b, LC_est_b);

% Plot estimated V_C and difference from true V_C
figure
plot(t, VC_b, 'r', t, V_C_est_b', 'b');
title('Comparison between V_{Cb} and V_{Cbest}');
xlabel('Time (s)');
ylabel('Magnitude');
legend('True output V_{Cb}', 'Estimated output V_{Cbest}');

figure
plot(t, VC_b - V_C_est_b, 'r');
title('Error in estimation of V_{Cb}');
xlabel('Time (s)');
ylabel('Magnitude');
legend('Difference');

% Plot estimated V_R and difference from true V_R
figure
plot(t, VR_b, 'r', t, V_R_est_b, 'b');
title('Comparison between V_{Rb} and V_{Rbest}');
xlabel('Time (s)');
ylabel('Magnitude');
legend('True output V_{Rb}', 'Estimated output V_{Rbest}');

figure
plot(t, VR_b - V_R_est_b, 'r');
title('Error in estimation of V_{Rb}');
xlabel('Time (s)');
ylabel('Magnitude');
legend('Difference');
