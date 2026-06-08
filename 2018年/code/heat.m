function u = heat(tmax, n, each, rho, C, k, d, u0, ue, k_e, k_s)
%HEAT 此处显示有关此函数的摘要
%   此处显示详细说明
dt = tmax/n;
dx = d./each;
m = sum(each);
u = zeros(n+1,m+1);                       % 行表示时间，列表示位置
u(1,:) = u0;                              % 初始零时刻全为u0,第i行表示(i-1)dt,j列表示(j-1)dx

% 构造矩阵
a = zeros(m,1); b = zeros(m+1,1); c = zeros(m,1); d = zeros(m+1,1);

% 设定边值条件
mu_e = (k_e / k(1)) * dx(1);
mu_s = (k_s / k(4)) * dx(4);
b(1) = 1 + mu_e;  c(1) = -1;  d(1) = mu_e * ue;
a(m) = -1;  b(m+1) = 1 + mu_s;  d(m+1) = mu_s * u0; 
for i = 2:n+1
    start = 2;                            % 除去边界,b从2开始到n
    finish = 0;                           
    for j = 1:4
        finish = finish + each(j);
        lamda = (k(j)/(C(j)*rho(j))) * (dt/(dx(j).^2));
        a(start-1:finish-1,1) = -lamda;
        b(start:finish,1) = 1 + 2*lamda;
        c(start:finish,1) = -lamda;
        d(start:finish,1) = u(i-1,start:finish);
        if j ~= 4
            a(finish) = - k(j)/dx(j);
            b(finish+1) = k(j)/dx(j) + k(j+1)/dx(j+1);
            c(finish+1) = - k(j+1)/dx(j+1);
            d(finish+1, 1) = 0;
        end
        start = start + each(j);
    end
    x = thomas(a,b,c,d);
    u(i,1:m+1) = x;
end
end

