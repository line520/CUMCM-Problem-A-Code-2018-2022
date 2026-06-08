clc, clear, close all;

% 拟合参数
k = [3.2570e-10; 4.9146e-11; 6.3837e-11; 7.9114e-11; 5.2674e-11; 2.5183e-11];
beta = [ 4.6366e-07; 2.5408e-06; 7.8538e-06; 1.5462e-06; 5.0409e-07; 2.6733e-06];

% 温区尺寸参数
L0 = 25; L1 = 30.5; l = 5;                          % 炉前炉后长度、每个小温区和间隙,单位cm
L = L0 + 11 * L1+ 10*l + L0;                        % 传送带总长度,单位cm

% 原器件尺寸参数
D = 0.15;                                           % 焊接区域厚度,单位mm

% 温度和速度初始值
T0 = 25; T5 = 25;                                   % T0表示环境温度,单位oC
v = 70/60;

rng(0); % 方便可以复现
% 遍历寻找最优曲线
dt = 0.05; dD = 0.005/1000;                          % 时间步长(s)、厚度步长(m)
m = round(D/(dD*1000));                              % 节点数

N1 = 5500; N2 = 6000; dN = 500; 
e = 0.05;                                            % 两次计算的误差
v1 = 0; v2 = 0; flag = true;
omega = 1e-8;                                        % 避免计算结果太好为零

while flag
    result = iter(T0, T5, L0, L1, l, D, k, beta, dt, dD, N1);
    m1 = max(result(:,6)); n1 = min(result(:,6));
    m2 = max(result(:,7)); n2 = min(result(:,7)); 
    [v1, ~] = min((result(:,6) - n1 + omega)/(m1 - n1) + ...
        (result(:,7) - n2 + omega)/(m2 - n2));
    fprintf('总个数为：%d\n', N1);
    fprintf('函数值为：%d\n', v1);

    result = iter(T0, T5, L0, L1, l, D, k, beta, dt, dD, N2);
    m1 = max(result(:,6)); n1 = min(result(:,6));
    m2 = max(result(:,7)); n2 = min(result(:,7)); 
    [v2, idx] = min((result(:,6) - n1 + omega)/(m1 - n1) + ...
        (result(:,7) - n2 + omega)/(m2 - n2));
    best_condition = result(idx,1:5);
    fprintf('总个数为：%d\n', N2);
    fprintf('函数值为：%d\n', v2);

    if abs(v1 - v2)/v2 < e 
        flag = false;
        fprintf('最佳面积为：%.6f\n', result(idx, 6));
        fprintf('最佳偏差为：%.6f\n', result(idx, 7));
        fprintf('最佳条件为：\n');
        fprintf('%.6f\n', best_condition);
        fprintf('指标值为：\n');
        fprintf('%d\n', v2);
    else
        fprintf('精度不足!加大次数！\n');
        N1 = N2;
        N2 = N2 + dN;
    end

end

function result = iter(T0, T5, L0, L1, l, D, k, beta, dt, dD, N)
    L = L0 + 11 * L1+ 10*l + L0;                        % 传送带总长度,单位cm
    m = round(D/(dD*1000));                             % 节点数

    randMatrix = rand(N,5);
    randMatrix(:,1) = 20*randMatrix(:,1) + 165;         % 小温区1~5
    randMatrix(:,2) = 20*randMatrix(:,2) + 185;         % 小温区6
    randMatrix(:,3) = 20*randMatrix(:,3) + 225;         % 小温区7
    randMatrix(:,4) = 20*randMatrix(:,4) + 245;         % 小温区8~9
    randMatrix(:,5) = 35*randMatrix(:,5) + 65;

    result = ones(N,7)*inf;
    po = 0;

    for i = 1:N
        T1 = randMatrix(i,1); T2 = randMatrix(i,2); T3 = randMatrix(i,3); 
        T4 = randMatrix(i,4); v = randMatrix(i,5)/60;
        t_total = L/v;
        T = [T0; T1; T2; T3; T4; T5; T0];
        u_t = heatmodel(T, L0, L1, l, v, D, k, beta, dt, dD);
        u = u_t(1 :(floor(t_total/dt)+1), ceil(m/2));
        if check(u, dt)                                 % 通过
            po = po + 1;
            result(po, 1:5) = randMatrix(i,:);
            S = get_square(u, dt);
            result(po, 6) = S;
            result(po, 7) = symmetry(u);
        end
    end
    fprintf('找到%d个解\n',po)
    result = result(1:po, :);
end

function square = get_square(u, dt)
    id1 = find(u >= 217, 1);
    [~, id2] = max(u);                  % 最大值索引
    
    bound_up = sum((u(id1:id2,1) - 217) * dt);
    bound_low = sum((u(id1-1:id2-1,1) - 217) * dt);

    if (bound_up - bound_low) > 2
        fprintf('精度不足，差距为：%.5f\n',bound_up - bound_low);
    end
    square = (bound_up + bound_low)/2;

end

function val = symmetry(u)
    id1 = find(u >= 217, 1);
    [~, id2] = max(u);                  % 最大值索引
    id3 = 2*id2 - id1;
    
    u1 = u(id1:1:id2);
    u2 = u(id3:-1:id2);

    val = sum((u1 - u2).^2,1); 
end

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