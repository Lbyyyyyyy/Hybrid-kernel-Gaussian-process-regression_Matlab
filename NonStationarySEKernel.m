% 非平稳线性平方指数核
function [K, Knm] = NonStationarySEKernel(XN, XM, theta)
    % 从 theta 中提取参数
    sigma_f = exp(theta(1)); % 信号方差（σ_f^2）
    a = theta(2);            % 线性调制函数的斜率（a）
    b = theta(3);            % 线性调制函数的截距（b）

    % 对输入应用线性调制函数，计算局部长度尺度
    l_XN = a*XN + b;
    l_XM = a*XM + b;

    % 计算SE核成分的平方距离
    sqdist_SE = pdist2(XN, XM).^2;

    % 计算调制后的长度尺度的平方
    l_sq = l_XN .* l_XM';
    
    % 非平稳SE核
    K = sigma_f^2 * exp(-sqdist_SE ./ (2 * l_sq));
    
    if nargout > 1
        % 如果需要，计算梯度
        
        % 关于 sigma_f 的梯度
        dK_dsigma_f = 2 * sigma_f * exp(-sqdist_SE ./ (2 * l_sq));

        % 关于 a 和 b 的梯度需要更复杂的处理，涉及到 l_XN 和 l_XM 的偏导数
        % 这里不提供 a 和 b 的梯度计算，因为它们需要考虑 l_XN 和 l_XM 关于 a 和 b 的偏导数

        % 暂时只返回 sigma_f 的梯度
        Knm = {[], dK_dsigma_f}; % 注意这里用空数组占位，表示没有计算关于a和b的梯度
    end
end
