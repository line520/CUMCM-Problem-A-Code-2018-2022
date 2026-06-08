function P = calP(rho)
    P_rho = [993160.69930514, -3262197.61451880, 4024401.29625003, -2207818.47823074, 454076.76110785];
    P = polyval(P_rho,rho);
    % if rho > 0.85
    %     drho = 0.0005;
    % else 
    %     drho = -0.0005;
    % end
    % 
    % P = 100;                                % 密度为0.85的时候压强为100MPa
    % rhonow = 0.85;
    % E_P = [0.00010004, -0.00108248, 5.47444434, 1531.86840585]; % 模拟出的E(P)的多项式系数
    % 
    % while abs(rhonow - rho) > abs(drho)
    %     E = polyval(E_P,P);
    %     dP = E/rho * drho;
    %     P = P + dP;
    %     rhonow = rhonow + drho;
    % end
end

