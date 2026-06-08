% 这个文件解决第一问
%% 编写说明
% 大框架：优化确定参数：误差函数->优化迭代
% heat利用差分形式求解各个时间各个地方的温度 输入：网格，常数，边界条件 输出：温度
% err函数表示误差
%% 求解代码
close all; clear,clc;
% 读入0s到5400s的假人外侧温度
data = readmatrix("CUMCM-2018-Problem-A-Chinese-Appendix.xlsx","Sheet","附件2","Range",'B3:B5403');
% 定义一些常数
rho = [300.00; 862.00; 74.20; 1.18];    % 1,2,3,4层密度(kg/m^3)
c = [1377; 2100; 1726; 1005];           % 1,2,3,4层比热容(J/kg*oC)
k = [0.082; 0.370; 0.045; 0.028];       % 1,2,3,4层热传导率(W/m*oC)
d = [0.6; 6; 3.6; 5]./1000;             % 1,2,3,4层厚度(mm) 2,4层为未知待定

% 划分网格数
tmax = 5400; 
dt = 0.5;
n = tmax/dt;                            % 划分时间段数
each = [60;600;360;500];                % 划分各层段数
% 边值条件
u0 = 37;                                % 初始各层的温度
ue = 75;                                % 初始环境温度
k_e0 = 50;                              % 初始参数


fun = @(k_e0) errfun(tmax, n, each, rho, c, k, d, u0, ue, k_e0, data);
%options = odeset('RelTol',1e-8,'AbsTol',1e-10);
k_e = fminbnd(fun,20,200);         % 求解出k_e和k_s
k_s = 1./((48.08-u0)/(ue-48.08) * (1./k_e + sum(d./k)));
fprintf('环境热交换系数为：%d\n 皮肤与4层热交换系数为：%d\n',k_e,k_s);
sol = heat(tmax, n, each, rho, c, k, d, u0, ue, k_e, k_s);

xmesh = 0:1e-5:sum(d);
tspan = 0: dt :tmax;
m = sum(each);
figure(1)
surf(xmesh,tspan,sol);
xlabel('x');
ylabel('t');
zlabel('u');
shading interp;

figure(2)
plot(xmesh, sol(n+1,1:m+1));
xlabel('x');
ylabel('u');

figure(3)
temp1 = sol(1:n+1,1);
temp2 = sol(1:n+1,each(1)+1);
temp3 = sol(1:n+1,each(1)+each(2)+1);
temp4 = sol(1:n+1,each(1)+each(2)+each(3)+1);
temp5 = sol(1:n+1,end);
plot(tspan,[temp1,temp2,temp3,temp4,temp5]);


function err = errfun(tmax, n, each, rho, c, k, d, u0, ue, k_e0, data)
k_s0 = 1./((48.08-u0)/(ue-48.08) * (1./k_e0 + sum(d./k)));
u = heat(tmax, n, each, rho, c, k, d, u0, ue, k_e0, k_s0);
m = sum(each);
tmp = u(1:n+1,m+1);
dt = tmax/n;
step = 1/dt; start = 1;
tmp_pre = zeros(length(data),1);
for i = 1:length(data)
    tmp_pre(i) = tmp(start);
    start = start + step;
end
err = sum((tmp_pre-data).^2)/2;
end
