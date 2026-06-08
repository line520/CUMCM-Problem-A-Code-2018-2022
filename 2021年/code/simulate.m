function R = simulate(p, r)
%SIMULATE 用来拟合最佳的旋转矩阵
% R 表示旋转矩阵
% r 表示基准圆半径
    theta0 = ones(3,1)*pi/4;
    fun = @(theta) err(p, theta, r);
    options = optimset('TolX',1e-10, 'MaxFunEvals', 1000, 'MaxIter', 5000);
    [theta_b, fval] = fminsearch(fun,theta0,options);
    R = Matrix(theta_b);
end

function val = err(p, theta, r)
    R = Matrix(theta);
    target = [0; 0; r];
    val = sqrt(sum((R*p - target).^2));
end

function R = Matrix(theta)
    Rx = [1     0      0;
          0     cos(theta(1))  sin(theta(1)); 
          0     -sin(theta(1)) cos(theta(1))];

    Ry = [cos(theta(2))     0     -sin(theta(2));
          0      1      0; 
          sin(theta(2))     0     cos(theta(2))];
    
    Rz = [cos(theta(3))     sin(theta(3))     0;
         -sin(theta(3))     cos(theta(3))     0; 
          0     0     1];
    R = Rx * Ry * Rz;
end

