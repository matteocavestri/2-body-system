function [M, C, K, h, x_g, A, V, rho, Fy, Fg] = inizialisation(c_pto, k_pto)
%INIZIALISATION c_pto [N/m/s]
%   Function to innizialize the problem
%   geometrical definitions and 

disp("Inizialisation started")

% Data geometry
% booth geometries are cylinders
global g rho_w rho_air w_w
syms t

%%%
% DEFINITION OF THE GEOMETRIC PARAMETER

r_1 = 1; % [m]
r_2up = 0.3; % [m]
h_1= 2;

A_1 = pi*(r_1^2-r_2up^2);
V_1 = A_1*h_1;

c_f = 1; % coefficent of dencity frequency
phi_0 = 0/4*(pi/c_f);

r_2dw(t) = 0.8 + (0.5*exp(1.i*c_f*(w_w*t + phi_0))); % [m] sinusoidal wave
% r_2dw(t) = 1.5 + 0.5*(sin(c_f*(w_w*t + phi_0)) + sin(c_f*(1/3*w_w*t + phi_0))/3 + sin(c_f*(1/5*w_w*t + phi_0))/5 + sin(c_f*(1/7*w_w*t + phi_0))/7); % [m] step wave

figure(1)
fplot(real(r_2dw)), grid on
title('moving radius')
xlabel('time [s]')
ylabel('radius [m]')
legend('r_{2-dw}')


figure("Name",'Rho')
subplot(211)
title('Mass analysis')
%fplot(r_2dw), grid on, hold on

h_2up = 3; % [m]
h_2dw = 1; % [m]

wide_st = 0.03; % 5cm wide
%%%

A_2up = pi*r_2up^2;

A_2dw = pi.*r_2dw.^2; % definition using sinusoidal functions

% Possibility of define A_2dw or r_2dw as a square wave function
V_2_up = A_2up*h_2up;
V_2_dw = A_2dw*h_2dw;
V_2 = V_2_dw + V_2_up;


h = [h_1, h_2up, h_2dw];
A = [A_1, A_2up, A_2dw];
V = [V_1, V_2];

% calculate the max weigth of the two bodies --> search maximum volume and
% calculate m_st using the max volume
rho_st = 2700; % steal density
t_0 = (-phi_0)/w_w;
V_st = [(h_1+2*wide_st)*pi*((r_1+wide_st)^2-(r_2up+wide_st)^2) - V_1, pi*((h_2dw+2*wide_st)*((r_2dw+wide_st)^2) + (h_2up+wide_st)*((r_2up+wide_st)^2)) - V_2]; % V_2(t_0) = max(V_2) 
m_st = rho_st.*V_st;

m_st = m_st(t_0);
m_sel = m_st*[0; 1];


% Data Dynamic System
a_12 = 0;
a_21 = 0;
m_1 = 3/2*m_st*[1; 0] + V_1*rho_air;
m_2 = m_st*[0; 1] + V_2*rho_air;

% x_g is a value containing the barycenter. It's neded to  the position an depth of the two body. 
x_g = [h_1/2; ((h_2up + h_2dw)/2*((h_2up+wide_st)*((r_2up+wide_st)^2)*rho_st + V_2_up*(rho_air-rho_st)) + (h_2dw/2*(pi*((h_2dw+2*wide_st)*((r_2dw+wide_st)^2))*rho_st + V_2_dw*(rho_air-rho_st))))/m_2];

fplot(real(m_1)), grid on, hold on
fplot(real(m_2)), grid on, hold on
ylabel('mass [kg]')
xlabel('time [s]')
legend("m_{1}", "m_{2}")

rho = [m_1/V_1, m_2/V_2];
rho_2_med = int(rho*[0; 1], [-pi pi])/(2*pi);

subplot(212)
title("Rho")
fplot(real(rho)), grid on, hold on,
fplot(real(rho_2_med)), grid on, hold off
ylabel('dencity [m^3/kg]')
xlabel('time [s]')
legend("rho_1", "rho_2", "rho_{2-med}")


M = [[m_1, a_12]; [a_21, m_2]]


% c_d = 1.15; % coefficent to unalize the drag dorce on cilinders 
c_1 = diff(rho_air*V_1, t); %{2*rho_w*A(1); % dumping factor related to the sharestress value}
c_2 = diff(rho_air*V_2, t); % component related to changing of body_2's volume  {%2*rho_w*A(2); % incert the value of the Drag force}


C = [[c_pto + c_1, -c_pto];[-c_pto, c_pto + c_2]]


% to define the rho matrix changing over time

k_1 = rho_w*g*A_1;
k_2 = rho_w*g*A_2up;

K = [[k_pto - k_1, -k_pto];[-k_pto, k_pto - k_2]]

Fy = rho_w.*g.*[A_1; A_2up]
Fg = -[m_1; m_2 - h_2dw.*rho_w*(A_2dw-A_2up)]

disp("Inizialisation ended")
end

