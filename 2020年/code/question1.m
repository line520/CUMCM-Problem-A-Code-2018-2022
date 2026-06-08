clear, clc, close all;

% 其余初始参数
% 各个温区参数设置
T0 = 25; T1 = 173; T2 = 198; T3 = 230; T4 = 257; T5 = 25;  % T0表示环境温度,单位oC  
T = [T0; T1; T2; T3; T4; T5; T0];
L0 = 25; L1 = 30.5; l = 5;                          % 炉前炉后长度、每个小温区和间隙,单位cm
v = 78/60;                                          % 传送带的速度,单位cm/s
D = 0.15;                                           % 焊接区域厚度,单位mm
k = [3.2570e-10; 4.9146e-11; 6.3837e-11; 7.9114e-11; 5.2674e-11; 2.5183e-11];
beta = [ 4.6366e-07; 2.5408e-06; 7.8538e-06; 1.5462e-06; 5.0409e-07; 2.6733e-06];

dt = 0.05; dD = 0.005/1000;                         % 时间步长(s)、厚度步长(m)
m = round(D/(dD*1000)); 
u = heatmodel(T, L0, L1, l, v, D, k, beta, dt, dD);

% 小温区3、6、7中心温度和小温区8结束处温度
fprintf('温区3中心温度：%.6f\n', u(floor((L0+2.5*L1+2*l)/(v*dt)), ceil(m/2)));
fprintf('温区6中心温度：%.6f\n', u(floor((L0+5.5*L1+5*l)/(v*dt)), ceil(m/2)));
fprintf('温区7中心温度：%.6f\n', u(floor((L0+6.5*L1+6*l)/(v*dt)), ceil(m/2)));
fprintf('温区8结束温度：%.6f\n', u(floor((L0+8*L1+7*l)/(v*dt)), ceil(m/2)));

% 绘制炉温曲线
L = L0 + 11 * L1+ 10*l + L0;                % 传送带总长度,单位cm
t_total = L/v;                              % 花费的总时间
pred = u(1: 0.5/dt :(floor(t_total/dt)+1), ceil(m/2));
time = 0:0.5:floor(t_total/dt)*dt;

figure(1)
plot(time', pred);

% 创建表格（第一行留空）
dataTable = array2table([time', pred],'VariableNames',{'时间(s)', '温度(摄氏度)'});

% 写入CSV文件（跳过第一行）
writetable(dataTable, 'result.csv', 'WriteMode', 'overwrite');
