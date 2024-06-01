%====绘制SOC的观测值与预测值散点图====
function print_OvsP(Ob,Pre,lim)
%++++ Input ++++
% Ob: m x 1 土壤有机碳的观测值  
% Pre: m x 1 土壤有机碳的预测值 
% lim: 表示x轴和Y轴的刻度上限
%
%++++ Output ++++
% 散点图
%-------------------------------------------------------------------------
% 创建散点图
figure
scatter(Ob, Pre,[],'k');
box on
grid on
grid minor
% 保持图像，以便在其上绘制线性回归线
hold on;

% 拟合线性回归模型
linear_fit = polyfit(Ob, Pre, 1);

% 计算拟合线的y值
fit_values = polyval(linear_fit, Ob);

% 绘制拟合线
plot(Ob, fit_values, 'r-', 'LineWidth', 2);

% 添加拟合线的方程和r值
% 计算相关系数和线性拟合的系数
r = corrcoef(Ob, Pre);
r = r(1,2);
str_eq = sprintf('y = %.4fx + %.4f', linear_fit(1), linear_fit(2));
str_R2 = sprintf('r = %.4f', r);
text(min(Ob), max(Pre), str_eq, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
text(min(Ob), max(Pre) - 5, str_R2, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% 绘制x=y参考线
hline = refline(1, 0);
set(hline, 'Color', 'k');

% 设置坐标轴标签
xlabel('SOC Observed (g/kg)');
ylabel('SOC Predicted (g/kg)');

% 设置坐标抽范围
xlim([0 lim]);
ylim([0 lim]);

% 设置图表标题
%title('Prediction vs Ground Truth of Cu Concentration');

% 取消保持状态
hold off;

end
