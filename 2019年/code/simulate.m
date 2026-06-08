clc, clear all;
data_1 = readmatrix("附件1-凸轮边缘曲线.xlsx","Range",'A2:B629');
% 这个是凸轮的极角和极径
data_2 = readmatrix('附件2-针阀运动曲线.xlsx','Range','A2:B46');
data_3 = readmatrix('附件2-针阀运动曲线.xlsx','Range','D2:E46');
% 这个是针阀运动的时间和高度
data_4 = readmatrix('附件3-弹性模量与压力.xlsx','Range','A2:B402');
% 这个是弹性模量和压力的数据

%% 首先对压力和弹性模量进行模拟
% 初始数据
p = data_4(:,1);
E = data_4(:,2);

% 最佳拟合
n = 4;
l = polyfit(p,E,n);
pred = polyval(l,p);
disp(num2str(l,'%.8f\t'));
err_last = sum((E-pred).^2);
disp(err_last);

figure(1)
plot(p,E); hold on;
plot(p,pred);
title('原始数据和拟合数据');

% 最佳拟合的由来
err = zeros(1,3);
for i = 3:5
    l = polyfit(p,E,i);
    pred = polyval(l,p);
    err(i-2) = sum((E-pred).^2);
end

figure(2)
plot(3:5,err);
title('不同次数多项式拟合效果');

%% 第二问的凸轮模拟
theta = data_1(:,1);
r = data_1(:,2);
[x,y] = pol2cart(theta, r);

figure(3)
polarplot(theta, r);
title('凸轮的形状');

figure(4)
plot(x,y); hold on;
title('凸轮');

% 凸轮拟合的思路：看图像类似于圆形，在基础上加上一个cos来进行拟合
var = cos(theta);
fun = @(a, v)  a(1) + a(2)*v;
[k,Q] = lsqcurvefit(fun, [0,0], var, r);
fprintf('%.6f\n', k);
[x_n, y_n] = pol2cart(theta, k(1) + k(2)*var);
plot(x_n,y_n);



%% 第二问针阀上升高度绘制
x1 = data_2(:,1); x2 = 0.45:0.01:2; x3 = data_3(:,1); % 对应的时刻
y1 = data_2(:,2); y2 = ones(size(x2)) * 2; y3 = data_3(:,2); % 对应针阀高度
x = [x1',x2,x3']; y = [y1',y2,y3'];
figure(5)
plot(x,y);
title('针阀上升高度');

%% 计算出P和rho之间的关系
% E_P = [0.00000035, -0.00003813, 0.01666712, 4.68786653, 1539.64578098];
rho = zeros(2001,1);
P = 0:0.1:200;
for i = 1:2001
    rho(i) = calrho_iter(P(i));
end

k = polyfit(rho,P,3);
disp(num2str(k,'%.8f\t'));
pred_P = polyval(k,rho);
pred_new = polyval([1.3316e5, -3.2602e5, 2.6814e5, -74048], rho);

q = polyfit(P,rho,2);
disp(num2str(q,'%.10f\t'));
pred_rho = polyval(q,P);

figure(6)
plot(rho,P); hold on;
plot(rho,pred_P);

figure(7)
plot(rho,P); hold on;
plot(rho, pred_new);

figure(8)
plot(P,rho); hold on;
plot(P,pred_rho);
% 计算得到最优参数为：
% P_rho = [993160.69930514, -3262197.61451880, 4024401.29625003, -2207818.47823074, 454076.76110785]
% rho_P = [-0.0000006582, 0.0005229512, 0.8042886771]

