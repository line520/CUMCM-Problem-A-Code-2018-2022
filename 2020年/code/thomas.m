function x = thomas(a, b, c, d)
%THOMAS 利用追赶法求解三对角代数方程
% 输入:
%   a: 次对角线元素（从第2行到第n行），长度 n-1
%   b: 主对角线元素，长度 n
%   c: 超对角线元素（从第1行到第n-1行），长度 n-1
%   d: 方程右侧向量，长度 n
% 输出:
%   x: 解向量

n = length(b);
x = zeros(n,1);
if length(a) ~= n-1 || length(c) ~= n-1 || length(d) ~= n
    error('输入向量维数不对！请重新输入！');
end

% 下面是追的过程，也就是化成上三角矩阵
for i = 2 : n
    m = a(i-1) / b(i-1);                % 用上一行消下一行的第一个元素需乘的倍数
    b(i) = b(i) - m * c(i-1);
    d(i) = d(i) - m * d(i-1);
end

% 下面是赶的过程，也就是根据上三角矩阵求解方程
x(n) = d(n) / b(n);
for i = n-1:-1:1
    x(i) = (d(i) - c(i) * x(i+1)) / b(i);
end

end

