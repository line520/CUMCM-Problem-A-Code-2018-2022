% 该文档主要进行数据预处理
% 生成数据如下：
%   relation 有关联的主索节点的索引数据
%   team     组成一个基准反射面的主索节点的索引数据
%   data1    基准球面上主索节点坐标
%   data2    促动器下端点(地锚点)坐标
%   data3    促动器基准态上端点(顶端)坐标

clc,clear, close all;
% 读取数据并转换为字符串数组
name = readtable('附件1.csv', 'Range', 'A2:A2227','ReadVariableNames',false);
name = string(name{:,:});
group = readtable('附件3.csv', 'Range', 'A2:C4301','ReadVariableNames',false);
group = string(group{:,:});
m = size(group, 1);
relation = zeros(3*m, 2);
team = zeros(m, 3);

for i = 1:m
    a = find(name == group(i, 1));
    b = find(name == group(i, 2));
    c = find(name == group(i, 3));
    
    relation(3*i-2, :) = [a, b];
    relation(3*i-1, :) = [a, c];
    relation(3*i, :) = [b, c];
    
    team(i, :) = [a, b, c];
end

relation = unique(relation, 'rows', 'sorted');
save('./data/name.mat','name');
save('./data/relation.mat','relation');
save('./data/team.mat','team');

% 读取坐标信息，按照name的顺序排序
d1 = readmatrix("附件1.csv",'Range','B2:D2227');
d2 = readmatrix("附件2.csv",'Range','B2:G2227');
data1 = d1(:,:);
data2 = d2(:,1:3);
data3 = d2(:,4:6);
save('./data/data1.mat','data1');
save('./data/data2.mat','data2');
save('./data/data3.mat','data3');
