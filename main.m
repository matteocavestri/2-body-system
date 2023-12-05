clear all
close all
global g rho_w rho_air w_w

% Data Environment
g = 9.81; % [m/s^2]
rho_w = 1000; % [kg/m^3]
rho_air = 1.2;

    % wave data
f_w = 0.3;
w_w = 2*pi*f_w;
amp_w = 1;

c_pto = 1000;
k_pto = 8000;

[M, C, K, h, x_g, A, V, rho, Fy, Fg] = inizialisation(c_pto, k_pto); % inizialise the model starting from geometry

% Create a function to optimize c_pto and k_pto over time.
%
% 

[y, x, x_d, x_dd] = solver(w_w, amp_w, M, C, K, Fy, Fg);

H_eval(w_w, M, C, K, Fy);

P = power_balancing(y, x, x_d, x_dd, A, h, x_g, c_pto);

fprintf("\nMechanism evaluated\n\nLoading plots\n")