clear, clc, close all;
% 读入某次实验数据
data = readmatrix('附件.xlsx','Range','A2:B710');
T_e = data(:,2);                                    % e表示实验数据                          

% 其余初始参数
% 各个温区参数设置
T0 = 25; T1 = 175; T2 = 195; T3 = 235; T4 = 255; T5 = 25;  % T0表示环境温度,单位oC  
T = [T0; T1; T2; T3; T4; T5; T0];
L0 = 25; L1 = 30.5; l = 5;                          % 炉前炉后长度、每个小温区和间隙,单位cm
v = 70/60;                                          % 传送带的速度,单位cm/s
D = 0.15;                                           % 焊接区域厚度,单位mm
k0 = ones(6,1) * 9e-11;
beta0 = ones(6,1) * (1e-6);
% k0 = ones(6,1);
% beta0 = ones(6,1) ;
x0 = [k0; beta0];

fun = @(x) errfun(T_e, T, L0, L1, l, v, D, x(1:6), x(7:12));
options = optimset('TolX',1e-8, 'PlotFcns',@optimplotfval, 'MaxFunEvals', 5000, 'MaxIter', 5000);
[xb, fval] = fminsearch(fun,x0,options);
fprintf('%.15f\n',xb);
disp(fval);

dt = 0.01; dD = 0.005/1000;                 % 时间步长(s)、厚度步长(m)
m = round(D/(dD*1000)); 
u = heatmodel(T, L0, L1, l, v, D, xb(1:6), xb(7:12), dt, dD);
pred = u(19/dt + 1: 0.5/dt :(19/dt + 1 +708*(0.5/dt)), ceil(m/2));


figure(1)
plot(data(:,1), data(:,2));hold on;
plot(data(:,1), pred);

function err = errfun(T_e, T, L0, L1, l, v, D, k, beta_t)
    if min(k) < 0 || min(beta_t) < 0
        err = inf;
        return
    end

    dt = 0.05; dD = 0.005/1000;                 % 时间步长(s)、厚度步长(m)
    m = round(D/(dD*1000)); 

    u = heatmodel(T, L0, L1, l, v, D, k, beta_t, dt, dD);
    pred = u(19/dt + 1: 0.5/dt :(19/dt + 1 +708*(0.5/dt)), ceil(m/2));
    
    err = sum((pred - T_e).^2);
end

% data = [
%     3.25701e-10;
%     4.9146e-11;
%     6.3837e-11;
%     7.9114e-11;
%     5.2674e-11;
%     2.5183e-11;
%     4.63663046e-7;
%     2.540828862e-6;
%     7.853769845e-6;
%     1.546215869e-6;
%     5.04088974e-7;
%     2.673317460e-6
% ];
