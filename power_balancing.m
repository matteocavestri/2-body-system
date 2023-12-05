function [P] = power_balancing(y, x, x_d, x_dd, A, h, x_g, c_pto)
%POWER_BALANCING Summary of this function goes here
%   Detailed explanation goes here

global g rho_w
syms t

disp('Power balancing started')

P_ext = abs(1/2.*c_pto.*x_d.^2); % power extracted

h_0 = [0, 1]*x - h*[0; 1; 1/2] %- y; % - [0, 1]*x_g - y;
p = - g*rho_w*(h_0);

figure(7)
sp(1) = subplot(211);
fplot(real(h_0))
xlabel('time [s]')
ylabel('meter [m]')
legend('h_{0_2}')

sp(1) = subplot(212);
fplot(real(p))
xlabel('time [s]')
ylabel('Pressure [Pa]')
legend('p_2')


% E = p.*(A*[0; 0; 1]).*(h*[0; 0; 1]);
P_in = p.*diff(A*[0; 0; 1], t).*(h*[0; 0; 1]); % power to move the mechanism

P = P_ext + P_in;
P_med = vpaintegral(P, [-pi pi])/(2*pi);
P_in_med = vpaintegral(P_in, [-pi pi])/(2*pi);

figure(6)
%sp(1) = subplot(211);
fplot(real(P_ext)), grid on, hold on,
fplot(real(P_in)), grid on, hold on,
fplot(real(P)), grid on, hold on,
fplot(real(P_med)), grid on, hold on,
fplot(real(P_in_med)), grid on, hold off,
xlabel('time [s]')
ylabel('Power [W]')
legend('P_{ext}', 'P_{in}', 'P_{net}', 'P_{med}', 'P_{in-med}')



%sp(2) = subplot(212);
%xlabel('time [s]')
%ylabel('Energy [J]')
%plot(real(E)), grid on, hold on,

%fplot(real(E)), grid on, hold on,
%fplot(imag(P_in)), grid on, hold on,
%fplot(imag(P)), grid on, hold off
%legend('P_{ext}', 'P_{in}', 'P_{net}')

disp('Power balancing ended')

end

