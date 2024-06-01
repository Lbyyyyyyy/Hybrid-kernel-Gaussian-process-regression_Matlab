% 局部多项式平方指数核
function [K, Knm] = LocalPolynomialSEKernel(XN, XM, theta)
    % 从 theta 中提取参数
    sigma_f = exp(theta(1)); % 信号方差（σ_f^2）
    sigma_l = exp(theta(2)); % 特征长度尺度（σ_l）

    % 设定调制参数a和b
    a = ones(size(XN, 2), 1); % 假设每个特征的权重相同
    b = 0;                    % 假设没有偏置

    % 计算调制后的特征内积 T
    T = (XN * a) * (XM * a)' + b; % 矩阵 XN 和 XM 每一行是一个样本点，每一列是一个特征

    % 局部多项式SE核
    K = sigma_f^2 * exp(-pdist2(XN * a, XM * a).^2 / (2 * sigma_l^2)) .* (1 + T).^2;
    
    if nargout > 1
        % 如果需要，计算梯度
        
        % 关于 sigma_f 的梯度
        dK_dsigma_f = 2 * K / sigma_f;
        
        % 关于 sigma_l 的梯度
        % 首先计算距离矩阵的每个元素关于 sigma_l 的导数
        XN_scaled = XN * a;
        XM_scaled = XM * a;
        pdist_scaled = pdist2(XN_scaled, XM_scaled).^2;
        dK_dsigma_l = K .* pdist_scaled / (sigma_l^3);
        
        % 返回计算出的梯度
        Knm = {dK_dsigma_l, dK_dsigma_f}; 
    end
end