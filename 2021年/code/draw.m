clc, clear, close all;
% 读取主索节点的数据
name1 = readtable('附件1.csv','Range','A2:A2227');
data1 = readmatrix("附件1.csv",'Range','B2:D2227');
p_x = data1(:,1); p_y = data1(:,2); p_z = data1(:,3); % 主索节点的坐标
r = sqrt(p_x.^2 + p_y.^2 + p_z.^2);
fprintf('平均长度：%.6f\n', mean(r));
fprintf('x最大：%.6f\n', max(abs(p_x)));
fprintf('y最大：%.6f\n', mean(abs(p_y)));
fprintf('z最大：%.6f\n', mean(abs(p_z)));

% 读取促动器数据
name2 = readtable('附件2.csv','Range','A2:A2227');
data2 = readmatrix("附件2.csv",'Range','B2:G2227');
a_x = data2(:,1); a_y = data2(:,2); a_z = data2(:,3); % 促动器下端点(地锚点)坐标
b_x = data2(:,4); b_y = data2(:,5); b_z = data2(:,6); % 促动器基准态上端点(顶端)坐标

if isequal(name1, name2)
    disp('两个字符串组相同');
end

% 读取反射面板主索节点编号
group = readtable('附件3.csv','Range','A2:C4301');

figure(1);
scatter3(p_x, p_y, p_z);
title('主索节点分布');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

figure(2);
scatter3(a_x, a_y, a_z);
title('地锚点分布');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

figure(3);
scatter3(b_x, b_y, b_z);
title('促动器顶端分布');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

figure(4);
scatter3(p_x, p_y, p_z);hold on;
scatter3(a_x, a_y, a_z);hold on;
scatter3(b_x, b_y, b_z); 
xlabel('X');
ylabel('Y');
zlabel('Z');