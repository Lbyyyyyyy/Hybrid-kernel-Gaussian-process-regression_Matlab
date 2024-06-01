% ====评估模型性能====
function PI = ModelAssessment(Yfit,Y,n,order)
% Model Assessment 
%
%---input---
% Yfit : Predicted value
% Y : Observed value
% n : number of sapmles
% order: In order to calculate RMSE,
%        if Yift and Y are the original value, order = 0; 
%        if Yift and Y are the logarithmic value, order = 1.
%
%---output---
% PI:Performance Index, a strcuct data,
%   rsaqured : coefficient of determination,0-1.
%   RMSE : root mean squared error
%   rRMSE : relative root mean squared error
%   RPD : the ratio of performance to deviation
%   RPIQ : the ratio of performance to interquartile range
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4;order = 0;end

TSS = sum((Y - mean(Y)).^2);
RSS = sum((Y - Yfit).^2);
rsquared = 1 - RSS / TSS;
RMSE = sqrt(RSS / n);
rRMSE = 100 * RMSE / mean(Y);
SD = sqrt(var(Y,1));
RPD = SD / RMSE;

Q1 = quantile(Y,0.25);
Q3 = quantile(Y,0.75);
IQ = Q3 - Q1;
RPIQ = IQ / RMSE;

if order == 0
    RMSE = RMSE;
else
    RMSE = sqrt(sum((exp(Yfit)-exp(Y)).^2) / n);
end

% -----Output
PI.R2 = rsquared;
PI.RMSE = RMSE;
% PI.RMSE = sqrt(sum((exp(Yfit)-exp(Y)).^2) / n);
PI.rRMSE = rRMSE;
PI.RPD = RPD;
PI.RPIQ = RPIQ;
end
