function rho = calrho_iter(P)
% 根据现在的P值计算对应密度
    if P > 100
        dP = 0.01;
    else 
        dP = -0.01;
    end
    
    rho = 0.85;                             % 100 MPa时候是0.85mg/mm^3
    Pnow = 100;
    % E_P = [0.00010004, -0.00108248, 5.47444434, 1531.86840585]; % 模拟出的E(P)的多项式系数
    E_P = [0.00000035, -0.00003813, 0.01666712, 4.68786653, 1539.64578098];
    while abs(Pnow - P) > abs(dP)
        E = polyval(E_P,Pnow);
        drho = rho/E * dP;
        rho = rho + drho;
        Pnow = Pnow + dP;
    end
    E = polyval(E_P, P);
end

% 计算得到 160MPa的时候密度为0.87113mg/mm^3
% 计算得到 0.5MPa的时候密度为0.804560mg/mm^3


