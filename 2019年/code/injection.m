data_2 = readmatrix('附件2-针阀运动曲线.xlsx','Range','A2:B46');
data_3 = readmatrix('附件2-针阀运动曲线.xlsx','Range','D2:E46');
x1 = data_2(:,1); x2 = 0.45:0.01:2; x3 = data_3(:,1);           % 对应的时刻
y1 = data_2(:,2); y2 = ones(size(x2)) * 2; y3 = data_3(:,2);    % 对应针阀高度
t_H = [x1',x2,x3']; H = [y1',y2,y3'];

dt = 0.01; t_w = 50;
result = zeros(200,1);
tnow = 0;
for i = 1:20000 
    tnow = tnow + dt;
    mark = round((tnow - floor(tnow./100) * 100 - t_w) / dt) + 1; 
    if mark <= 33 && mark > 0
        A_p = pi * ( (1.25 + H(mark) * tan(pi/20)).^2 - 1.25.^2);
    elseif mark > 33 && mark <= 213
        A_p = pi * 0.7 * 0.7;
    elseif mark > 213 && mark <= 246
        A_p = pi * ( (1.25 + H(mark) * tan(pi/20)).^2 - 1.25.^2);
    else
        A_p = 0;
     end
    result(i) = A_p;
end
plot(0.01:0.01:200,result);
