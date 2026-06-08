clc, clear, close all;
%% 基本参数
% 浮子基本参数,后缀为f
m_f = 4866;                     % 质量，单位为kg
R_f = 1;                        % 底部半径，单位m
h_f1 = 3;                       % 圆柱部分高度,单位m
h_f2 = 0.8;                     % 圆锥部分高度，单位m

% 振子基本参数，后缀为v
m_v = 2433;                     % 质量，单位为kg

% 其他参数
w = 2.2143;                     % 入射波浪频率，单位为s-1
m_a = 1165.992;                 % 垂荡附加质量，单位kg
lamda_h = 167.8395;             % 垂荡兴波阻尼系数，单位N*s/m
f_d = 4890;                     % 垂荡激励力振幅，单位N
k1 = 80000;                     % 弹簧刚度，单位N/m
rho = 1025;                     % 海水的密度，单位kg/m^3
g = 9.8;                        % 重力加速度，单位m/s^2
h_s = 2.00;                     % 由第一问得到的静止时水线高度，单位m

% 下面采用变步长搜索法
lamda_k_max = 100000;
alpha_max = 1;
T = 2 * pi / w; 
t_total = 600*T;                % 总时间
tspan = linspace(0,t_total,10000);
dt = t_total/10000;             % 时间步长
y0 = [0, 0, 0, 0];              % 初始值

% 首先大步长搜索
lamda_k = 0;                    % 初始值
dlamda = 500;                   % lamda_1 步长
dalpha = 0.1;                   % alpha   步长
n = lamda_k_max/dlamda;
m = alpha_max/dalpha;
lamda_k_best = inf;
alpha_best = inf;
val_max = 0;
result = zeros(n+1,m+1);

for i = 1:n+1
    alpha = 0;
    for j = 1:m+1
        fun = @(t, y) odefun(y, t, rho, g, R_f, h_f2, lamda_k, lamda_h, m_f, m_a, m_v, k1, f_d, w, h_s, alpha);
        [t, y] = ode45(fun, tspan, y0);
        lamda_1 = lamda_k*(abs(y(:,2) - y(:,4))).^alpha;
        val = (sum(lamda_1.*(y(:,2) - y(:,4)).^2 * dt))/ t_total;
        result(i,j) = val;
        if val > val_max
            val_max = val;
            lamda_k_best = lamda_k;
            alpha_best = alpha;
        end
        alpha = alpha + dalpha;
    end
    lamda_k = lamda_k + dlamda;
end

fprintf('最佳阻尼系数为：%.6f N*s/m\n',lamda_k_best);
fprintf('对应功率为：%.6f W\n', val_max);
fprintf('对应幂为：%.6f\n', alpha_max);

figure(6)
surf(0:dalpha:alpha_max, 0:dlamda:lamda_k_max, result);

% 接下来小步长搜索
lamda_k = lamda_k_best - dlamda;% 初始值
alpha = alpha_best - dalpha;
dlam = 0.5 ;                       % 步长
dalp = 0.01;
n = 2*dlamda/dlam;  
m = 2*dalpha/dalp;

for i = 1:n+1
    for j = 1:m+1
        fun = @(t, y) odefun(y, t, rho, g, R_f, h_f2, lamda_k, lamda_h, m_f, m_a, m_v, k1, f_d, w, h_s, alpha);
        [t, y] = ode45(fun, tspan, y0);
        val = (sum(abs(lamda_k*(y(:,2) - y(:,4)).^(2+alpha) * dt)))/ t_total;
        if val > val_max
            val_max = val;
            lamda_k_best = lamda_k;
            alpha_best = alpha;
        end
        alpha = alpha + dalp;
    end
    lamda_k = lamda_k + dlam;
end

fprintf('最佳阻尼系数为：%.6f N*s/m\n',lamda_k_best);
fprintf('对应功率为：%.6f W\n', val_max);
fprintf('对应幂为：%.6f\n', alpha_max);

fun = @(t, y) odefun(y, t, rho, g, R_f, h_f2, lamda_k_best, lamda_h, m_f, m_a, m_v, k1, f_d, w, h_s, alpha_best);
[t, y] = ode45(fun, tspan, y0);
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

function dydt = odefun(y, t, rho, g, R_f, h_f2, lamda_k, lamda_h, m_f, m_a, m_v, k1, f_d, w, h_s, alpha)
% y 存储了各个因变量按顺序是
% 浮子：x_f, v_f 振子：x_v, v_v
    if y(1) >= -h_s
        k_s = rho * g * pi * R_f^2;
    else
        k_s = (((h_f2 + y(1) + 2)/ h_f2) * R_f).^2 * pi;
    end
    lamda_1 = lamda_k*(abs(y(2) - y(4))^alpha);
    dydt = zeros(4,1);
    dydt(1) = y(2);
    dydt(2) = (-lamda_h*y(2) + lamda_1*(y(4)-y(2)) - k_s*y(1) ...
        + k1*(y(3) - y(1)) + f_d*cos(w*t))/(m_f + m_a);
    dydt(3) = y(4);
    dydt(4) = - lamda_1/m_v * (y(4) - y(2)) - k1/m_v * (y(3) - y(1));
end
