clc, clear, close all;
% 这里的直线阻尼器系数不为常数
%% 基本参数
% 浮子基本参数,后缀为f
m_f = 4866;                     % 质量，单位为kg
R_f = 1;                        % 底部半径，单位m
h_f1 = 3;                       % 圆柱部分高度,单位m
h_f2 = 0.8;                     % 圆锥部分高度，单位m

% 振子基本参数，后缀为v
m_v = 2433;                     % 质量，单位为kg

% 其他参数
w = 1.4005;                     % 入射波浪频率，单位为s-1
m_a = 1335.535;                 % 垂荡附加质量，单位kg
lamda_h = 656.3616;             % 垂荡兴波阻尼系数，单位N*s/m
f_d = 6250;                     % 垂荡激励力振幅，单位N
k1 = 80000;                     % 弹簧刚度，单位N/m
rho = 1025;                     % 海水的密度，单位kg/m^3
g = 9.8;                        % 重力加速度，单位m/s^2

% 首先计算静止时候，水线位置
V_z = 1/3 * pi * R_f^2 * h_f2;  % 圆锥部分体积
V_p = (m_f + m_v) / rho;        % 排水体积
h_s = (V_p - V_z)/(pi * R_f * R_f); % 静止时水线位置

lamda_k = 10000;                % 直线阻尼器的比例系数

fun = @(t, y) odefun(y, t, rho, g, R_f, h_f2, lamda_k, lamda_h, m_f, m_a, m_v, k1, f_d, w, h_s);
tspan = 0:0.2:180;
y0 = [0, 0, 0, 0];
[t, y] = ode45(fun, tspan, y0);

writematrix(y,'result1-2.xlsx','Range','B3');
writematrix(tspan','result1-2.xlsx','Range','A3');
figure(1)
plot(t, y(:,1));
xlabel('t');
ylabel('浮子位移');

figure(2)
plot(t, y(:,2));
xlabel('t');
ylabel('浮子速度');

figure(3)
plot(t, y(:,3));
xlabel('t');
ylabel('振子位移');

figure(4)
plot(t, y(:,4));
xlabel('t');
ylabel('振子速度');
function dydt = odefun(y, t, rho, g, R_f, h_f2, lamda_k, lamda_h, m_f, m_a, m_v, k1, f_d, w, h_s)
% y 存储了各个因变量按顺序是
% 浮子：x_f, v_f  振子：x_v, v_v 
    if y(1) >= -h_s
        k_s = rho * g * pi * R_f^2;
    else
        k_s = (((h_f2 + y(1) + 2)/ h_f2) * R_f).^2 * pi;
    end
    lamda_1 = lamda_k*(abs(y(2) - y(4))^0.5);
    dydt = zeros(4,1);
    dydt(1) = y(2);
    dydt(2) = (-lamda_h*y(2) + lamda_1*(y(4)-y(2)) - k_s*y(1) ...
        + k1*(y(3) - y(1)) + f_d*cos(w*t))/(m_f + m_a);
    dydt(3) = y(4);
    dydt(4) = - lamda_1/m_v * (y(4) - y(2)) - k1/m_v * (y(3) - y(1));
end



