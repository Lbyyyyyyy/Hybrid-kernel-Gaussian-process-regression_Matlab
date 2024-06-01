%====ͨ����������ϵ��ȷ��PLSRǱ��������====
function LV = SelectNumberofLV(X,Y,LV_max)
%++++ Input ++++
% X: m x n  (Sample matrix)
% Y: m x 1  (measured property) ԭʼֵ
% LV_max: the max number of LVs
%
%++++ Output ++++
% LV:a strcuct data
% 
[m,n] = size(X);
[~,~,~,~,~,pctVar,mse] = plsregress(X,log(Y),LV_max,'CV',10);
figure; % ���Ʋ�ͬ���ɷ���������������Ľ��ͳ̶�
plot(1:LV_max,cumsum(100*pctVar(2,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in Y');
figure; % ���Ʋ�ͬ���ɷ������ľ��������
plot(0:LV_max,mse(2,:),'b-o');
xlabel('Number of components');
ylabel('Estimated Mean Squared Prediction Error');

ADrsquared = [];
TSS = sum((log(Y) - mean(log(Y))).^2);
for i = 1:LV_max
    [~,~,~,~,beta] = plsregress(X,log(Y),i,'CV',10); % ��log(Y)�������
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
LV.ADrsquared = ADrsquared; % ������R2
LV.ADrsquared_max = ADrsquared_max; % ���������R2
LV.optLV = AD_index; % ��������R2��Ӧ��LV����
LV.note2 = '*** The following is based on global min MSE + 1SD'
LV.RMSE = RMSE; % ���������
LV.RMSE_min = RMSE_min % ��ͬLV����СRMSE
LV.RMSE_min_1SD = RMSE(indexSD); % ����СRMSE��һ����׼���ڵ�...
                                 ...��СLV����Ӧ��RMSE
LV.ADrsquared_max_1SD = ADrsquared(indexSD); % ����СRMSE��һ����׼���ڵ�...
                                             ...��СLV����Ӧ�ĵ�����R2
                                                 
LV.opt_1SD = indexSD;
end