clear, clc, close all;

% 各个温区参数设置
T0 = 25; T1 = 182; T2 = 203; T3 = 237; T4 = 254; T5 = 25;  % T0表示环境温度,单位oC 
T = [T0; T1; T2; T3; T4; T5; T0];
L0 = 25; L1 = 30.5; l = 5;                          % 炉前炉后长度、每个小温区和间隙,单位cm
L = L0 + 11 * L1+ 10*l + L0;                        % 传送带总长度,单位cm
D = 0.15;                                           % 焊接区域厚度,单位mm
dt = 0.1; dD = 0.001/1000;                          % 时间步长(s)、厚度步长(m)
m = round(D/(dD*1000));                             % 节点数

% 拟合参数
k = [3.2570e-10; 4.9146e-11; 6.3837e-11; 7.9114e-11; 5.2674e-11; 2.5183e-11];
beta = [ 4.6366e-07; 2.5408e-06; 7.8538e-06; 1.5462e-06; 5.0409e-07; 2.6733e-06];

% 遍历求解最大速度
v0 = 70/60; vmax = 80/60;                         % 最大和最小传送带的速度,单位cm/s
dv = 0.01/60;                                         % 步长,单位cm/s
n = (vmax - v0)/dv;
vb = 0;
for i = 1:n+1
    t_total = L/v0;                               % 花费的总时间
    u_t = heatmodel(T, L0, L1, l, v0, D, k, beta, dt, dD);
    u = u_t(1: 1 :(floor(t_total/dt)+1), ceil(m/2));
    if check(u, dt) && vb < v0   % 通过
        vb = v0;
    end
    v0 = v0 + dv;
end

fprintf('最佳的运行速率为：%.5fcm/min', vb*60);

% 绘制炉温曲线
u_t = heatmodel(T, L0, L1, l, vb, D, k, beta, dt, dD);
t_total = L/vb;                              % 花费的总时间
u = u_t(1: 0.5/dt :(floor(t_total/dt)+1), ceil(m/2));
time = 0:0.5:floor(t_total/dt)*dt;

figure(1)
plot(time', u);

function flag = check(u, dt)
% 这个函数用来检查是否满足制程界限
    flag = false;
    u = u(u>=30);
    d = diff(u);
    if (max(u) >= 240 && max(u) <= 250) && (max(abs(d)) <= 3 * dt)
        % 判断上升下降斜率 和 峰值温度
        u1 = u(u>217);
        t1 = size(u1,1) * dt;                % 大于217度时间

        d1 = [d; -1];
        u2 = u(d1 >= 0);                   % 表示在上升阶段
        u3 = u2(u2 >= 150 & u2 <= 190);
        t3 = size(u3,1) * dt;
        if t1 >= 40 && t1 <= 90 && t3 >= 60 && t3 <= 120
            flag = true;
        end
    end
end