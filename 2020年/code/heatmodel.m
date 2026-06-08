function u = heatmodel(T, L0, L1, l, v, D, k, beta_t, dt, dD)
% T0表示环境温度
% T 表示各个小温区的温度
% k 表示各个区域的传热系数
% beta_t 表示各个区域的换热系数

    L = L0 + 11 * L1+ 10*l + L0;                % 传送带总长度,单位cm
    t_total = L/v;                              % 花费的总时间
    dx = v * dt/100;                            % 距离步长(m)
    n = ceil(t_total/dt);                       % 时间点个数
    m = round(D/(dD*1000));                     % 距离点个数
    u = zeros(n+1, m+1);                        % 行表示时间，列表示距离

    % 初始条件
    u(1,:) = T(1);
    xnow = 0;                                   % 从零位置开始出发
    
    % 构建方程组的向量
    a = zeros(m,1); b = zeros(m+1,1); c = zeros(m,1); d = zeros(m+1,1);
    eta = 5;
    % flag = false;

    for i = 2:n+1
        xnow = xnow + dx*100;                       % 传送带运行到的位置

        [T_s, alpha, beta] = getparameter(xnow, T, k, beta_t, L0, L1, l, dt, dD, eta);

        a(m) = -beta/(dD+beta); a(1:m-1, 1) = -alpha;
        b(1) = 1; b(m+1) = 1; b(2:m,1) = 2*alpha + 1;
        c(1) = -beta/(dD+beta); c(2:m,1) = -alpha;
        d(1) = dD/(dD+beta) * T_s; d(m+1) = dD/(dD+beta) * T_s;
        d(2:m) = u(i-1,2:m);

        u(i,:) = thomas(a, b, c, d);

        
        % if abs(max(u(i,ceil(m/2))- u(i-1,ceil(m/2)))) > 3 * dt
        %     flag = true;
        %     break
        % end

    end
    % if flag 
    %     pred = ones(709,1) * inf;
    % else
    %    pred = u(19/dt + 1: 0.5/dt :(19/dt + 1 +708*(0.5/dt)), ceil(m/2));
    % end
end

