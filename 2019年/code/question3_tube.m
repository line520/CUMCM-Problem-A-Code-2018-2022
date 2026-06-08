clc, clear, close all;
% 初始数据
% 这个是针阀运动的时间和高度
data_2 = readmatrix('附件2-针阀运动曲线.xlsx','Range','A2:B46');
data_3 = readmatrix('附件2-针阀运动曲线.xlsx','Range','D2:E46');
x1 = data_2(:,1); x2 = 0.45:0.01:2; x3 = data_3(:,1);           % 对应的时刻
y1 = data_2(:,2); y2 = ones(size(x2)) * 2; y3 = data_3(:,2);    % 对应针阀高度
t_H = [x1',x2,x3']; H = [y1',y2,y3'];

C = 0.85;
A = pi * 0.7 * 0.7;                         % 输油管的横截面积,单位mm^2
B = pi * 0.7 * 0.7;                         % 减压阀的横截面积,单位mm^2
S = pi * 2.5 * 2.5;                         % 柱塞的横截面积，单位mm^2
V = pi * 5 * 5 * 500;                       % 高压油管的体积,单位mm^3
P_0 = 100;                                  % 高压油管初始压强,单位：MPa

d = 1e-5;                                   % 每次搜索的步长
dt = 5;
dw = 7.4e-4;                                % 2.5e-2 rad/ms 
n = 5;                                     % 搜索次数
tmax = 10000; 
valmin = inf;   dwmin = inf; twmin = inf;
for i = 1:n+1
    t_w = 35;
    for j = 1:n+1
        P = fun(P_0, S, C, A, V, dw, H, tmax, t_w, B);
        val = sum(abs(P - 100));
        disp(val);
        if val < valmin
            valmin = val;
            dwmin = dw;
            twmin = t_w;
        end
        t_w = t_w + dt;
    end
    dw = dw + d;
end

P = fun(P_0, S, C, A, V, dwmin, H, 10000, twmin, B);
dt = 0.01;
figure(1)
plot(0:dt:10000,P);
title('最佳结果');
fprintf('搜索结果：\n 最佳角速度：%.6f \n %.6f\n目标函数:%.6f', dwmin*100, twmin, valmin);
%% 自定义函数区
function P = fun(P_0, S, C, A, V, dw, H, tmax, t_w, B) 
    % 内置参数
    dt = 0.01;                             % 时间步长
    tnow = 0;                              % 当前时刻
    wnow = 3.14;                           % 假设从最低点开始
    n = tmax/dt;

    P = zeros(1,n+1);
    P(1) = P_0;                            % 初始时刻是P_0
    rho_l = 0.804560;                      % 计算出来0.5MPa时的密度
    m_c = 20*rho_l + (7.239 - 2.413) * S * rho_l; % 初始时刻柱塞内充满低压燃油
    % m_c = 0;
    % 构建凸轮点列 r = 4.8260 + 2.4130 * cos(theta)
    theta = linspace(0,2*pi,4000);
    r = 4.8260 + 2.4130 * cos(theta);
    x = r .* cos(theta); y = r .* sin(theta);
    flag = true;                               % 让底部加油只加一次

    for i = 2:n+1
        tnow = tnow + dt;
        wnow = mod(wnow + dw, 2*pi); 
        % 计算流入的油量
        h = max(x*cos(wnow) - y*sin(wnow));     % 柱塞升程
        
        if abs(h - 2.413) < 1e-3 && flag        % 降到最低处补充燃油
            m_c = 20*rho_l + (7.239 - 2.413) * S * rho_l;
            flag = false;
        elseif h > 3 && ~flag
            flag = true;
        end
        V_c = 20 + (7.239 - h) * S;
        rho_e = m_c / V_c;
        
        % 如果左边压力比右边大则单向阀开启
        P_c = calP(rho_e);
        if P_c > P(i-1) && rho_e > 0.804288
            R = C * A * sqrt(2 * (P_c - P(i-1)) * rho_e) * dt;
            m_c = m_c - R;
        else 
            R = 0;
        end

        % 计算喷出油的质量

        % 第一个喷油嘴
        [E,rho] = calrho(P(i-1));
        % 计算喷油面积
        mark1 = round((mod(tnow, 100)) / dt) + 1;
        if mark1 <= 33 && mark1 > 0
            A_p1 = pi * ( (1.25 + H(mark1) * tan(pi/20)).^2 - 1.25.^2);
        elseif mark1 > 33 && mark1 <= 213
            A_p1 = pi * 0.7 * 0.7;
        elseif mark1 > 213 && mark1 <= 246
            A_p1 = pi * ( (1.25 + H(mark1) * tan(pi/20)).^2 - 1.25.^2);
        else
            A_p1 = 0;
        end

        if P(i-1) > 0
            Q1 = C * A_p1 * sqrt(2 * P(i-1) * rho) * dt; 
        else 
            Q1 = 0;
        end
        
        % 第二个喷油嘴 延迟了t_w秒进行喷油
        mark2 = round((mod(tnow, 100) - t_w) / dt) + 1;
        if mark2 <= 33 && mark2 > 0
            A_p2 = pi * ( (1.25 + H(mark2) * tan(pi/20)).^2 - 1.25.^2);
        elseif mark2 > 33 && mark2 <= 213
            A_p2 = pi * 0.7 * 0.7;
        elseif mark2 > 213 && mark2 <= 246
            A_p2 = pi * ( (1.25 + H(mark2) * tan(pi/20)).^2 - 1.25.^2);
        else
            A_p2 = 0;
        end

        if P(i-1) > 0
            Q2 = C * A_p2 * sqrt(2 * P(i-1) * rho) * dt; 
        else 
            Q2 = 0;
        end

        % 减压阀
        % 控制策略：当压强超过一定值的时候开启减压阀
        if P(i-1) > 102
            Q3 = C * B * sqrt(2 * P(i-1) * rho) * dt;
        else
            Q3 = 0;
        end

        P(i) = P(i-1) + E/(rho * V) * (R - Q1 - Q2 - Q3);


        if P(i) < 0
            P(i) = 0;
        end
    end
end
