%====通过调整决定系数确定PLSR潜变量数量====
function LV = SelectNumberofLV(X,Y,LV_max)
%++++ Input ++++
% X: m x n  (Sample matrix)
% Y: m x 1  (measured property) 原始值
% LV_max: the max number of LVs
%
%++++ Output ++++
% LV:a strcuct data
% 
[m,n] = size(X);
[~,~,~,~,~,pctVar,mse] = plsregress(X,log(Y),LV_max,'CV',10);
figure; % 绘制不同主成分数量对于因变量的解释程度
plot(1:LV_max,cumsum(100*pctVar(2,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in Y');
figure; % 绘制不同主成分数量的均方根误差
plot(0:LV_max,mse(2,:),'b-o');
xlabel('Number of components');
ylabel('Estimated Mean Squared Prediction Error');

ADrsquared = [];
TSS = sum((log(Y) - mean(log(Y))).^2);
for i = 1:LV_max
    [~,~,~,~,beta] = plsregress(X,log(Y),i,'CV',10); % 对log(Y)进行拟合
    yfit_logOC = [ones(m,1) X]*beta;
    RSS = sum((log(Y) - yfit_logOC).^2);
    ADrsquared(i) = 1-((m-1)/(m - LV_max -1))*(RSS/TSS);
    RMSE(i) = sqrt(sum((exp(yfit_logOC)-Y).^2) / m);
end
RMSE_min = min(RMSE);
[ADrsquared_max,AD_index] = max(ADrsquared);
indexSD = find(RMSE<=RMSE_min + std(RMSE));
indexSD = min(indexSD);
% ++++Output
LV.note1 = '*** The following is based on global max ADrsquared'
LV.ADrsquared = ADrsquared; % 调整后R2
LV.ADrsquared_max = ADrsquared_max; % 调整后最大R2
LV.optLV = AD_index; % 最大调整后R2对应的LV数量
LV.note2 = '*** The following is based on global min MSE + 1SD'
LV.RMSE = RMSE; % 均方根误差
LV.RMSE_min = RMSE_min % 不同LV的最小RMSE
LV.RMSE_min_1SD = RMSE(indexSD); % 在最小RMSE的一个标准差内的...
                                 ...最小LV所对应的RMSE
LV.ADrsquared_max_1SD = ADrsquared(indexSD); % 在最小RMSE的一个标准差内的...
                                             ...最小LV所对应的调整后R2
                                                 
LV.opt_1SD = indexSD;
end