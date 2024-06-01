% 定义高斯过程回归的核函数，返回核矩阵和参数梯度
% 定义混合核函数，即将Squared Exponential Kernel和Matern32 Kernel相加
% 四个参数
function [K, Knm] = HybridKernel(XN, XM, theta)
    % Parameters for Squared Exponential Kernel
    sigma_f_se = exp(theta(2)); % Scaling factor for Squared Exponential Kernel
    sigma_l_se = exp(theta(1)); % Length scale for Squared Exponential Kernel
    
    % Parameters for Matern32 Kernel
    sigma_f_m = exp(theta(4)); % Scaling factor for Matern32 Kernel
    sigma_l_m = exp(theta(3)); % Length scale for Matern32 Kernel
    
    % Calculate squared distances for the Squared Exponential Kernel
    sqdist_se = pdist2(XN, XM).^2;
    % Calculate the Squared Exponential Kernel
    K_se = sigma_f_se^2 * exp(-0.5 * sqdist_se / sigma_l_se^2);
    
    % Calculate Euclidean distances for the Matern32 Kernel
    r_m = pdist2(XN, XM);
    % Calculate the Matern32 Kernel
    K_matern32 = sigma_f_m^2 * (1 + sqrt(3) * r_m / sigma_l_m) .* exp(-sqrt(3) * r_m / sigma_l_m);
    
    % Combine the kernels
    K = K_se + K_matern32;
    
    if nargout > 1
        % Calculate gradients for the Squared Exponential Kernel
        dK_dsigma_f_se = 2 * sigma_f_se * exp(-0.5 * sqdist_se / sigma_l_se^2);
        dK_dsigma_l_se = (sqdist_se .* K_se) / (sigma_l_se^3);
        
        % Calculate gradients for the Matern32 Kernel
        dK_dsigma_f_m = 2 * sigma_f_m * (1 + sqrt(3) * r_m / sigma_l_m) .* exp(-sqrt(3) * r_m / sigma_l_m);
        dK_dsigma_l_m = sigma_f_m^2 * (sqrt(3) * r_m .* exp(-sqrt(3) * r_m / sigma_l_m)) .* ((-3) * r_m / sigma_l_m^2 + 3 / sigma_l_m);
        
        % Combine the gradients
        Knm = {dK_dsigma_l_se + dK_dsigma_l_m, dK_dsigma_f_se + dK_dsigma_f_m};
    end
end


