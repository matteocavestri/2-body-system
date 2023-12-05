function [y, x, x_d, x_dd] = solver(w_w, amp_w, M, C, K, Fy, Fg)
syms t
disp("Solver started")
    % define the movement of the mechanism
Lambda = -w_w.^2.*M + 1i.*w_w.*C + K;

H = Lambda\Fy;

x_0 = Lambda\Fg;
x_01 = [1, 0]*x_0;
x_02 = [0, 1]*x_0;

y = amp_w.*exp(1i.*w_w*t);

x = H.*y + x_0;
x_12 = [1,-1]*x;

figure(3)

sp(1)=subplot(211);
fplot(imag(y)), grid on, hold on
fplot(real(x_01)), grid on, hold on
fplot(real(x_02)), grid on, hold off
legend('wave', 'x_0_1', 'x_0_2')
title('Sea wave function')
ylabel('y wave [m]')
xlabel('time [s]')

sp(2)=subplot(212);
fplot(real(x)), grid on, hold on
fplot(real(x_12)), grid on, hold off
legend('x_1', 'x_2', 'x_{1-2}')
title('Model responce')
ylabel('x model [m]')
xlabel('time [s]')


x_d = diff(x_12, t);
x_dd = diff(x_d, t);

figure(4)
title('Velocity and Aceleration')

sp(1) = subplot(211);
fplot(real(x_d)), grid on
legend('Velocity')
ylabel('v_{pto} [m/s]')
xlabel('time [s]')

sp(2)=subplot(212);
fplot(real(x_dd)), grid on
legend('Aceleration')
ylabel('a_{pto} [m/s^2]')
xlabel('time [s]')
disp("Solver ended")


