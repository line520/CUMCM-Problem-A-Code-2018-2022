function p = position(alpha, beta, R)
% 用来计算天体位置
    if alpha >= - pi/2 || alpha <= pi/2
        x = sqrt((R * cos(beta))^2/(1+(tan(alpha)^2)));
    else
        x = -sqrt((R * cos(beta))^2/(1+(tan(alpha)^2)));
    end

    if alpha >= 0 || alpha <= pi
        y = sqrt((R * cos(beta) * tan(alpha))^2/(1+(tan(alpha)^2)));
    else 
        y = -sqrt((R * cos(beta) * tan(alpha))^2/(1+(tan(alpha)^2)));
    end

    if beta >= 0 || beta <= pi
        z = sqrt((R * sin(beta))^2);
    else
        z = -sqrt((R * sin(beta))^2);
    end
    p = [x; y; z];
end

