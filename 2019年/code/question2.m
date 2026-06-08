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
S = pi * 2.5 * 2.5;                         % 柱塞的横截面积，单位mm^2
V = pi * 5 * 5 * 500;                       % 高压油管的体积,单位mm^3
P_0 = 100;                                  % 高压油管初始压强,单位：MPa

d = 1e-6;                                   % 每次搜索的步长
dt = 5;
dw = 2.7e-4;                                % 2.5e-2 rad/ms 
n = 5;                                     % 搜索次数
tmax = 10000; 
valmin = inf;   dwmin = inf; twmin = inf;
for i = 1:n+1
    t_w = 30;
    for j = 1:n+1
        P = fun(P_0, S, C, A, V, dw, H, tmax, t_w);
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

P = fun(P_0, S, C, A, V, dwmin, H, 10000, twmin);
dt = 0.01;
figure(1)
plot(0:dt:10000,P);
title('最佳结果');
fprintf('搜索结果：\n 最佳角速度：%.6f \n %.6f\n目标函数:%.6f', dwmin*100, twmin, valmin);
%% 自定义函数区
function P = fun(P_0, S, C, A, V, dw, H, tmax, t_w) 
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
    
    hope = true;
    record = 0;

    for i = 2:n+1
        tnow = tnow + dt;
        wnow = mod(wnow + dw, 2*pi); 
        % 计算流入的油量
        h = max(x*cos(wnow) - y*sin(wnow));     % 柱塞升程
        
        if abs(h - 2.413) < 1e-3 && flag        % 降到最低处补充燃油 abs(h - 2.413) < 1e-3
            %m_c = m_c + (7.239 - 2.413) * S * rho_l;
            m_c = 20*rho_l + (7.239 - 2.413) * S * rho_l;
            flag = false;
            % disp('use');
        elseif h > 3 && ~flag
            flag = true;
            % disp('reset');
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
        
        % if hope && abs(wnow - 2*pi) > dw
        %     record = record + R;
        % elseif abs(wnow - 2*pi) < dw
        %     record = record + R;
        %     hope = false;
        %     fprintf('旋转一周：%.6f\n', record);
        % end

        % 计算喷出油的质量
        % 计算喷油面积
        mark = round((mod(tnow, 100) - t_w) / dt) + 1;
        [E,rho] = calrho(P(i-1));
        if mark <= 33 && mark > 0
            A_p = pi * ( (1.25 + H(mark) * tan(pi/20)).^2 - 1.25.^2);
        elseif mark > 33 && mark <= 213
            A_p = pi * 0.7 * 0.7;
        elseif mark > 213 && mark <= 246
            A_p = pi * ( (1.25 + H(mark) * tan(pi/20)).^2 - 1.25.^2);
        else
            A_p = 0;
        end

        if P(i-1) > 0
            Q = C * A_p * sqrt(2 * P(i-1) * rho) * dt; 
        else 
            Q = 0;
        end
        
        if hope && abs(wnow - 2*pi) > dw
            record = record + Q;
        elseif abs(wnow - 2*pi) < dw
            record = record + Q;
            hope = false;
            fprintf('旋转一周：%.6f\n', record);
        end

        P(i) = P(i-1) + E/(rho * V) * (R - Q);


        if P(i) < 0
            P(i) = 0;
        end
    end
end
