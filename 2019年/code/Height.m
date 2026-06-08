clc, clear, close all;
% 导入数据
data_1 = readmatrix("附件1-凸轮边缘曲线.xlsx","Range",'A2:B629');
theta_h = data_1(:,1);   r_h = data_1(:,2);

% 内置参数
wnow = 0;
dw = 0.0001;
n = ceil(2*pi/dw);
result = zeros(n,1);

for i = 1:n+1
    % 得到凸轮的函数 r = 4.8260 + 2.4130 * cos(theta)
    theta = linspace(0,2*pi,1000);
    r = 4.8260 + 2.4130 * cos(theta);
    x = r .* cos(theta); y = r .* sin(theta);
    h = max(x*cos(wnow) - y*sin(wnow));
    result(i) = h;
    wnow = wnow + dw; 
end

fprintf('最大上升高度：%.6f',max(result));
fprintf('最小上升高度：%.6f',min(result));

p1 = plot(theta_h,r_h); hold on;
p2 = plot(0:dw:n*dw, result);
title('凸轮极径和柱塞上升高度对比');
legend([p1, p2], {'极径', '上升高度'});
