function x = thomas_optimized(a, b, c, d)
% THOMAS 利用追赶法求解三对角代数方程（优化版）
% 输入:
%   a: 次对角线元素（从第2行到第n行），长度 n-1
%   b: 主对角线元素，长度 n
%   c: 超对角线元素（从第1行到第n-1行），长度 n-1
%   d: 方程右侧向量，长度 n
% 输出:
%   x: 解向量

n = length(b);
% 输入维度检查
if length(a) ~= n-1 || length(c) ~= n-1 || length(d) ~= n
    error('输入向量维数不对！请重新输入！');
end

% 预处理：复制输入防止修改原始数据
b_mod = b(:);  % 确保列向量
d_mod = d(:);
c_mod = c(:);
a_mod = a(:);

% 设置数值稳定性容差
tol = 1e-12;

% ==== 追过程（前向消元）====
for i = 2 : n
    % 检查主元避免除零
    if abs(b_mod(i-1)) < tol
        error('主对角线元素 %.6f 接近零，算法不稳定!', b_mod(i-1));
    end
    
    % 计算乘子并更新当前行
    m = a_mod(i-1) / b_mod(i-1);
    b_mod(i) = b_mod(i) - m * c_mod(i-1);
    d_mod(i) = d_mod(i) - m * d_mod(i-1);
end

% ==== 赶过程（后向替换）====
x = zeros(n, 1);
% 最后一行单独处理
if abs(b_mod(n)) < tol
    error('消元后主对角线元素 %.6f 接近零！', b_mod(n));
end
x(n) = d_mod(n) / b_mod(n);

% 倒序求解
for i = n-1 : -1 : 1
    % 使用已存储的c_mod值
    x(i) = (d_mod(i) - c_mod(i) * x(i+1)) / b_mod(i);
end

end

