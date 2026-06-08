clc, clear, close all;
%% 导入数据
name = load('./data/name.mat').name;
relation = load('./data/relation.mat').relation;
data1 = load('./data/data1.mat').data1;                 % 基准球面上主索节点坐标
data2 = load('./data/data2.mat').data2;                 % 促动器下端点(地锚点)坐标
data3 = load('./data/data3.mat').data3;                 % 促动器基准态上端点(顶端)坐标

%% 其他参数
% 装置参数
R = 300.4;                                          % 圆面半径，单位m
F = 0.466 * R;                                      % 焦点位置，单位m
d = 150;                                            % 抛物面口径的一半，单位m
alpha = 36.795/180*pi; beta = 78.169/180*pi;        % 天体位置
u = 0.07/100;                                       % 节点距离变化幅度

%% 准备工作
% 理想抛物面
epsilon = 0.36588;
p = 2*(F + epsilon);   f = R + epsilon;

% 旋转矩阵
location = position(alpha, beta, R);                       % 此时天体所在位置
Rot_o = simulate(location, R);                             % 输出旋转矩阵
Rot = inv(Rot_o);                                   % 将主索节点转回原位置

%  旋转并筛选
p_a = data1 * Rot_o';                              % 旋转后的主索节点坐标
p_c = data3 * Rot_o';                              % 旋转后促动器上端坐标
x_lim = sqrt((R*d)^2/(d^2 + (d^2/(2*p) - f)^2));
condition = (p_a(:,1) >= -x_lim) & (p_a(:,1) <= x_lim);

n = sum(condition);                                % 筛选出的节点数量
name_use = name(condition);
index = zeros(n,1); j = 1;
for i = 1:size(condition)
    if condition(i)
        index(j) = i;
        j = j + 1;
    end
end

% 提前计算各个节点之间的距离来提高速度
adj_matrix = sparse(size(p_a, 1), size(p_a, 1));
for i = 1:size(relation, 1)
    id_i = relation(i, 1);
    id_j = relation(i, 2);
    dist = sum((p_a(id_i,:) - p_a(id_j,:)).^2);
    adj_matrix(id_i, id_j) = dist;
    adj_matrix(id_j, id_i) = dist;
end
adj_matrix = adj_matrix(index, index);

% 预计算连接关系（只计算上三角部分）
[row, col] = find(triu(adj_matrix > 0));
distance0 = full(adj_matrix(sub2ind([n, n], row, col))); % 存储距离平方

% 预计算D_val和l_val
D = sqrt(sum(p_c(index,:).^2, 2));             % 促动器上端点到原点的距离
l = sum((p_a(index,:) - p_c(index,:)).^2, 2);  % 主索节点到促动器上端点的距离平方
%% 优化过程
p_total0 = p_a(condition,:);                       % 优化主索节点初始值
delta0 = ones(n,1)*0.1;                               % 促动器的收缩量
x0 = [p_total0(:); delta0];                        % 初始值      

fun = @(x) objective(x, p, f, n);
options = optimoptions('fmincon', 'Display', 'iter', ...
                      'Algorithm', 'interior-point', ...
                      'MaxIterations', 500, ...    % 先减少迭代次数测试
                      'FunctionTolerance', 1e-6, ...
                      'StepTolerance', 1e-10, ...
                      'OptimalityTolerance', 1e-6, ...
                      'MaxFunctionEvaluations', 1000000, ...
                      'UseParallel', true);
% 约束条件
A = []; b = []; Aeq = []; beq = [];                 % 线性约束（空）
lb = [-300.5*ones(3*n,1); -0.6*ones(n,1)];             % 变量下界
ub = [300.5*ones(3*n,1); 0.6*ones(n,1)];               % 变量上界
nonlcon = @(x) Constraint(x, index, p_c, u, n, row, col, distance0, D, l); % 非线性约束函数

[x_opt, fval] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

% 重塑优化结果以便后续使用
p_opt = reshape(x_opt(1:3*n), n, 3);               % 优化后的主索节点坐标
delta_opt = x_opt(3*n+1:end);                      % 优化后的促动器收缩量

function err = objective(x, p, f, n)
    p_total = reshape(x(1:3*n), n, 3);
    x_p = p_total(:,1);
    y_p = p_total(:,2);
    z_p = p_total(:,3);

    err = sum(((x_p.^2 + y_p.^2)/(2*p) - f - z_p).^2);
end

function [c, ceq] = Constraint(x, index, p_c, u, n, row, col, distance0, D, l)
% c     是非线性不等式约束
% ceq   是非线性等式约束

% u     表示节点之间距离变化幅度

    p_total = reshape(x(1:3*n), n, 3);
    delta = x(3*n+1:end);
    
    % 计算当前距离平方（向量化）
    diff = p_total(row,:) - p_total(col,:);
    distance= sum(diff.^2, 2);
    
    % 不等式约束：距离变化不超过u%
    c = abs(distance - distance0) - (u^2) * distance0;
    
    % 等式约束：促动器长度约束
    ratio = (D + delta) ./ D;
    x_d = ratio .* p_c(index,:);
    ceq = sum((p_total - x_d).^2, 2) - l;
end


