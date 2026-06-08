function [T_s, alpha, beta] = getparameter(xnow, T, k, beta_t, L0, L1, l, dt, dD, eta)
% 获取此时环境参数
    if xnow < L0 - eta                           % 炉前
        T_s = T(1);
        alpha = k(1) * dt/(dD).^2;
        beta = beta_t(1);
    elseif xnow < L0 + eta
        ratio = (xnow - L0 + eta )/(2*eta);
        T_s = T(1) + ratio * (T(2) - T(1));
        alpha = (k(1) + ratio * (k(2) - k(1))) * dt/(dD).^2;
        beta = beta_t(1) + ratio * (beta_t(2) - beta_t(1));
    elseif xnow < 5*L1 + 4*l + L0 - eta           % 温区1~5
        T_s = T(2);
        alpha = k(2) * dt/(dD).^2;
        beta = beta_t(2);
    elseif xnow < 5*L1 + 5*l + L0 + eta
        ratio = (xnow - 5*L1 - 4*l - L0 + eta)/(l+2*eta);
        T_s = T(2) + ratio * (T(3) - T(2));
        alpha = (k(2) + ratio * (k(3) - k(2))) * dt/(dD).^2;
        beta = beta_t(2) + ratio * (beta_t(3) - beta_t(2));
    elseif xnow < 6*L1 + 5*l + L0 - eta           % 温区6
        T_s = T(3);
        alpha = k(3) * dt/(dD).^2;
        beta = beta_t(3);
    elseif xnow < 6*L1 + 6*l + L0 + eta
        ratio = (xnow - 6*L1 - 5*l - L0 + eta)/(l+2*eta);
        T_s = T(3) + ratio * (T(4) - T(3));
        alpha = (k(3) + ratio* (k(4) - k(3))) * dt/(dD).^2;
        beta = beta_t(3) + ratio * (beta_t(4) - beta_t(3));
    elseif xnow < 7*L1 + 6*l + L0 - eta         % 温区7
        T_s = T(4);
        alpha = k(4) * dt/(dD).^2;
        beta = beta_t(4);
    elseif xnow < 7*L1 + 7*l + L0 + eta
        ratio = (xnow - 7*L1 - 6*l - L0 + eta)/(l+2*eta);
        T_s = T(4) + ratio * (T(5) - T(4));
        alpha = (k(4) + ratio* (k(5) - k(4))) * dt/(dD).^2;
        beta = beta_t(4) + ratio * (beta_t(5) - beta_t(4));
    elseif xnow < 9*L1 + 8*l + L0 - eta          % 温区8~9
        T_s = T(5);
        alpha = k(5) * dt/(dD).^2;
        beta = beta_t(5);
    elseif xnow < 9*L1 + 9*l + L0 + eta
        ratio = (xnow - 9*L1 - 8*l - L0 + eta )/(l+2*eta);
        T_s = T(5) + ratio * (T(6) - T(5));
        alpha = (k(5) + ratio* (k(6) - k(5))) * dt/(dD).^2;
        beta = beta_t(5) + ratio * (beta_t(6) - beta_t(5));
    elseif xnow < 11*L1 + 10*l + L0 - eta      % 温区10~11
        T_s = T(6);
        alpha = k(6) * dt/(dD).^2;
        beta = beta_t(6);
    else                                    % 炉后
        T_s = T(1);
        alpha = k(6) * dt/(dD).^2;
        beta = beta_t(6);
    end
end

