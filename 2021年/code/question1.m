clc, clear, close all;
% 读取数据
name = readtable('附件1.csv','Range','A2:A2227');
data1 = readmatrix("附件1.csv",'Range','B2:D2227');
data2 = readmatrix("附件2.csv",'Range','B2:G2227');
group = readtable('附件3.csv','Range','A2:C4301');
p_x = data1(:,1); p_y = data1(:,2); p_z = data1(:,3); % 主索节点的坐标
a_x = data2(:,1); a_y = data2(:,2); a_z = data2(:,3); % 促动器下端点(地锚点)坐标
b_x = data2(:,4); b_y = data2(:,5); b_z = data2(:,6); % 促动器基准态上端点(顶端)坐标

% 圆面数据
R = 300.4;  F = 0.466 * R; eta = 0.6;

epsilon_min = 0; epsilon_max = 0.6;                % 偏移量的最大值和最小值
epsilon = 0; n = 5000;
de = (epsilon_max - epsilon_min)/n;
valmin = inf; epsmin = inf;
% x = p_x; z = p_z;
theta = linspace(pi + acos(150/R), 2*pi - acos(150/R),10000);
x = (R*cos(theta))';
z = (R*sin(theta))';

for i = 1:n
    val = fun(x, z, F, R, epsilon, eta);
    if val < valmin
        valmin = val;
        epsmin = epsilon;
    end
    epsilon = epsilon + de;
end

fprintf('最佳函数值%.6f\n',valmin);
fprintf('对应偏差为：%.5f\n', epsmin);

function val = fun(x, z, F, R, epsilon, eta)
    m = size(x,1);
    x_p = zeros(m, 1);
    z_p = zeros(m, 1);
    for i = 1:m
        if x(i) == 0
            continue;
        end
        k = z(i)/x(i);
        coff = [1, -4*k*(F+epsilon), -4*(F+epsilon)*(R + epsilon)];
        x_root = roots(coff);
        z_root = (x_root.^2)./(4*(F+epsilon)) - R - epsilon;
        for j = 1:2
            if abs(x(i)/z(i) - x_root(j)/z_root(j)) < 1e-6
                x_p(i) = x_root(j);
                z_p(i) = z_root(j);
            end
        end
    end
    if max( sqrt((x_p - x).^2 + (z_p -z).^2) ) > eta
        val = inf;
        return
    end
    val = sqrt(sum((x_p - x).^2 + (z_p -z).^2));
end

