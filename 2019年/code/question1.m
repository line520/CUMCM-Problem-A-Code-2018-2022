clc ,clear all;
% 初始数据
C = 0.85;
A = pi * 0.7 * 0.7;                         % 输油管的横截面积,单位mm^2
V = pi * 5 * 5 * 500;                       % 高压油管的体积,单位mm^3
P_0 = 100;                                  % 高压油管初始压强,单位：MPa
P_f = 160;                                  % 输油管阀右侧压强,单位：MPa

% 下面的代码是经历调试后得到的代码，在这之前经过反复调整的最佳时间总是在0.3左右
% 原因是和最佳图形之间不匹配，最佳的是喷油的在中间的，所以在下面的代码中，我们引
% 入了一个非常重要的假设：喷油在每个100ms中间，入油在开头就有，对应参数为t_w

t_w = 48.5;                                 % 表示开始喷油的时刻
t = 0.283;                                  % 假设初始的单向管的开启时间,单位ms
dt = 0.001;                                 % 每次搜索的步长
n = 10;                                      % 搜索次数
valmin = inf;   tmin = inf; twmin = inf;
for i = 1:n+1
    for j = 1:n+1
        P = fun(P_0, P_f, C, A, V, t, dt, t_w);
        val = sum(abs(P - 100));
        if val < valmin
            valmin = val;
            tmin = t;
            twmin = t_w;
        end
        t_w = t_w + 0.5;
    end
    t = t + dt;
end

% 搜索最佳的喷油时间
% t_w = 40;
% dt = 0.5;
% n = 20;
% twmin = inf; valmin = inf;
% for i = 1: n+1
%     P = fun(P_0, P_f, C, A, V, t, dt, t_w);
%     val = sum(abs(P - 100));
%     if val < valmin
%         valmin = val;
%         twmin = t_w;
%     end
%     t_w = t_w + dt;
% end


P = fun(P_0, P_f, C, A, V, tmin, dt, twmin);
% dt = 0.01;
figure(1)
plot(0:dt:1000,P);
title('最佳结果');
fprintf('搜索结果：\n 最佳时间：%.6f \n %.6f\n 目标函数:%.6f',tmin, twmin, valmin);
%% 自定义函数区
function P = fun(P_0, P_f, C, A, V, t, dt, t_w)  
    % dt = 0.01;
    t_o = t;                               % t是单向阀开启的时间
    t_s = 0;                               % 表示暂停的时间
    tmax = 1000;                           % 我们用300ms时间来看函数值大小
    tnow = 0;                              % 当前时刻
    n = tmax/dt;
    P = zeros(1,n+1);
    P(1) = P_0;                            % 初始时刻是P_0
    rho_t = 0.87113;                        % 计算出来160MPa时的密度

    for i = 2:n+1
        tnow = tnow + dt;
        % 计算流入的油量
        if t_s > 0 && t_s >= dt
            t_s = t_s - dt;
            R = 0;                              % 处于没有开启状态，流入为零
        elseif t_s > 0 && t_s < dt
            t_s = 0;
            t_o = t_o - (dt - t_s);
            R = C * A * sqrt(2 * (P_f - P(i-1)) / rho_t) * (dt - t_s) * rho_t; % 流入油量的质量
        elseif t_o > 0 && t_o >= dt
            R = C * A * sqrt(2 * (P_f - P(i-1)) / rho_t) * dt * rho_t; % 流入油量的质量
            t_o = t_o - dt;
        elseif t_o > 0 && t_o < dt 
            R = C * A * sqrt(2 * (P_f - P(i-1)) / rho_t) * t_o * rho_t; % 流入油量的质量
            t_s = 10 - (dt - t_o);
            t_o = t;
        else
            t_s = 10 - dt;
            t_o = t;
            R = 0;
        end
        

        % 计算喷出油的质量, 从t_w时刻开始喷油
        [E, rho] = calrho(P(i-1));
        mark = tnow - floor(tnow./100) * 100 - t_w;   % 利用mark来判断现在所处阶段

        if mark > 0 && mark <= 0.2
            Q = 100 * (mark - dt/2) * dt * rho;
        elseif mark > 0.2 && mark <= 2.2
            Q = 20 * dt * rho;
        elseif mark > 2.2 && mark <= 2.4
            Q = (20 - 100 * (mark - 2.2 - dt/2)) * dt * rho;
        else 
            Q = 0;
        end
        
        P(i) = P(i-1) + E/(rho * V) * (R - Q);
    end
end
