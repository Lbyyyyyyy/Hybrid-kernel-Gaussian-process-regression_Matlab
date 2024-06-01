% 定义高斯过程回归的核函数，返回核矩阵和参数梯度
% 定义混合核函数，即将Squared Exponential Kernel和Matern32 Kernel相加
% 两个参数
function [K, Knm] = HybridKernel_two(XN, XM, theta)
    % 两个核共享的参数
    sigma_f = exp(theta(2)); % 尺度因子，共享
    sigma_l = exp(theta(1)); % 长度尺度，共享
    
    % 计算平方指数核的平方距离
    sqdist_se = pdist2(XN, XM).^2;
    % 计算平方指数核
    K_se = sigma_f^2 * exp(-0.5 * sqdist_se / sigma_l^2);
    
    % 计算Matern32核的欧氏距离
    r_m = pdist2(XN, XM);
    % 计算Matern32核
    K_matern32 = sigma_f^2 * (1 + sqrt(3) * r_m / sigma_l) .* exp(-sqrt(3) * r_m / sigma_l);
    
    % 结合两个核
    K = K_se + K_matern32;
    
    if nargout > 1
        % 如有需要，计算共享参数的梯度
        dK_dsigma_f = 2 * sigma_f * (exp(-0.5 * sqdist_se / sigma_l^2) + (1 + sqrt(3) * r_m / sigma_l) .* exp(-sqrt(3) * r_m / sigma_l));
        dK_dsigma_l = sigma_f^2 * (-0.5 * sqdist_se ./ sigma_l^3 .* exp(-0.5 * sqdist_se / sigma_l^2) + (sqrt(3) * r_m .* exp(-sqrt(3) * r_m / sigma_l)) .* ((-3) * r_m / sigma_l^3 + 3 / sigma_l^2));
        
        % 结合梯度
        Knm = {dK_dsigma_l, dK_dsigma_f};
    end
end
%----------------------------------------------------------------
% 1.theta(1) 和 theta(2) 分别是两个核共享的长度尺度和尺度因子。
% 2.平方指数核（Squared Exponential Kernel）和 
% Matern 3/2核都使用相同的长度尺度 sigma_l 和尺度因子 sigma_f。
% 当需要梯度信息时，我们计算两个核对于这些共享参数的梯度，然后相加。
% 这个改写使得核函数更简洁，并且核参数更容易解释，
% 因为现在两个核使用的是相同的长度尺度和尺度因子。
% 这在某些情况下可能是合理的，尤其是当你假设不同的核捕捉到的空间相关性具有相似的特征时。
% 





