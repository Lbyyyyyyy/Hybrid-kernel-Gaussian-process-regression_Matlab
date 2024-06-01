% 不同核函数的高斯过程回归模型
%% 1.导入数据
X_data_input_Train = xlsread('LUCAS_Spectra_2015_Agriculture_Train.xlsx');
Y_data_input_Train = xlsread('LUCAS_Topsoil_2015_Agriculture_Train.xlsx');
X_data_input_Test = xlsread('LUCAS_Spectra_2015_Agriculture_Test.xlsx');
Y_data_input_Test = xlsread('LUCAS_Topsoil_2015_Agriculture_Test.xlsx');

% 训练数据
Spectra_Train = X_data_input_Train(2:end,54:end); % 500-2498nm
OC_Train = Y_data_input_Train(:,10);

% 测试数据
Spectra_Test = X_data_input_Test(2:end,54:end); % 500-2498nm
OC_Test = Y_data_input_Test(:,10);

%% 2.预处理
% 训练数据中，光谱数据的一阶导数
Spectra_Train_FD = diff(Spectra_Train,1,2) / 2;
% 测试数据中，光谱数据的一阶导数
Spectra_Test_FD = diff(Spectra_Test,1,2) / 2;

% 计算训练集和测试集中样本数量 和 光谱波段数量
[n_Train,p_Train] = size(Spectra_Train_FD);
[n_Test,p_Test] = size(Spectra_Test_FD);

%% 3.基于不同核函数的高斯过程回归模型
% 基于matern32核-------------------------------------------------------
rng default
Mdl_gpr_matern32 = fitrgp(Spectra_Train_FD,log(OC_Train),...
    'BasisFunction','constant',...
    'KernelFunction','matern32',...
    'OptimizeHyperparameters','auto', ...
    'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus',...
    'Kfold',10))
tic
yfit_Train_gpr_matern32 = predict(Mdl_gpr_matern32,Spectra_Train_FD);
yfit_Test_gpr_matern32 = predict(Mdl_gpr_matern32,Spectra_Test_FD);
toc
% 性能评估
PI_Train_gpr_matern32 = ModelAssessment...
    (yfit_Train_gpr_matern32,log(OC_Train),n_Train,1);
PI_Test_gpr_matern32 = ModelAssessment...
    (yfit_Test_gpr_matern32,log(OC_Test),n_Test,1);

% 基于squared核------------------------------------------------------------
rng default
Mdl_gpr_squared = fitrgp(Spectra_Train_FD,log(OC_Train),...
    'BasisFunction','constant',...
    'KernelFunction','squaredexponential',...
    'OptimizeHyperparameters','auto', ...
    'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus',...
    'Kfold',10,...
    'MaxObjectiveEvaluations',10))
tic
yfit_Train_gpr_squared = predict(Mdl_gpr_squared,Spectra_Train_FD);
yfit_Test_gpr_squared = predict(Mdl_gpr_squared,Spectra_Test_FD);
toc
% 性能评估
PI_Train_gpr_squared = ModelAssessment...
    (yfit_Train_gpr_squared,log(OC_Train),n_Train,1);
PI_Test_gpr_squared = ModelAssessment...
    (yfit_Test_gpr_squared,log(OC_Test),n_Test,1);
% 基于sqared_1 ---------------------------------------------------
load('Mdl_gpr_squared_1.mat')
tic
yfit_Train_gpr_squared_1 = predict(Mdl_gpr_squared_1,Spectra_Train_FD);
yfit_Test_gpr_squared_1 = predict(Mdl_gpr_squared_1,Spectra_Test_FD);
toc
% 性能评估
PI_Train_gpr_squared_1 = ModelAssessment...
    (yfit_Train_gpr_squared_1,log(OC_Train),n_Train,1);
PI_Test_gpr_squared_1 = ModelAssessment...
    (yfit_Test_gpr_squared_1,log(OC_Test),n_Test,1);

% 基于ex核---------------------------------------------------------
rng default
Mdl_gpr_ex = fitrgp(Spectra_Train_FD,log(OC_Train),...
    'BasisFunction','constant',...
    'KernelFunction','exponential',...
    'OptimizeHyperparameters','auto', ...
    'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus',...
    'Kfold',10))
tic
yfit_Train_gpr_ex = predict(Mdl_gpr_ex,Spectra_Train_FD);
yfit_Test_gpr_ex = predict(Mdl_gpr_ex,Spectra_Test_FD);
toc
% 性能评估
PI_Train_gpr_ex = ModelAssessment...
    (yfit_Train_gpr_ex,log(OC_Train),n_Train,1);
PI_Test_gpr_ex = ModelAssessment...
    (yfit_Test_gpr_ex,log(OC_Test),n_Test,1);

% 基于rational核-----------------------------------------------------------
rng default
Mdl_gpr_rational = fitrgp(Spectra_Train_FD,log(OC_Train),...
    'BasisFunction','constant',...
    'KernelFunction','rationalquadratic',...
    'OptimizeHyperparameters','auto', ...
    'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus',...
    'Kfold',10))
tic
yfit_Train_gpr_rational = predict(Mdl_gpr_rational,Spectra_Train_FD);
yfit_Test_gpr_rational = predict(Mdl_gpr_rational,Spectra_Test_FD);
toc
% 性能评估
PI_Train_gpr_rational = ModelAssessment...
    (yfit_Train_gpr_rational,log(OC_Train),n_Train,1);
PI_Test_gpr_rational = ModelAssessment...
    (yfit_Test_gpr_rational,log(OC_Test),n_Test,1);

%% 4. 混合核函数的高斯过程回归
% 四个参数的形式
Hybrid = @(X,Y,theta)HybridKernel(X,Y,theta);
theta0=[1.5, 0.2, 1.5, 0.2];

rng default
Mdl_gpr_Hybrid = fitrgp(Spectra_Train_FD,log(OC_Train),...
    'BasisFunction','constant',...
    'KernelFunction',Hybrid,...
    'KernelParameters',theta0,...
    'OptimizeHyperparameters','auto', ...
    'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus',...
    'Kfold',10))
tic
yfit_Train_gpr_Hybrid = predict(Mdl_gpr_Hybrid,Spectra_Train_FD);
yfit_Test_gpr_Hybrid = predict(Mdl_gpr_Hybrid,Spectra_Test_FD);
toc
% 性能评估
PI_Train_gpr_Hybrid = ModelAssessment...
    (yfit_Train_gpr_Hybrid,log(OC_Train),n_Train,1);
PI_Test_gpr_Hybrid = ModelAssessment...
    (yfit_Test_gpr_Hybrid,log(OC_Test),n_Test,1);

%---------------------------------------------------------------
% 两个参数的形式
Hybrid_two = @(X,Y,theta)HybridKernel_two(X,Y,theta);
theta0=[mean(std(Spectra_Train_FD)), std(log(OC_Train))/sqrt(2)];

rng default
Mdl_gpr_Hybrid_two = fitrgp(Spectra_Train_FD,log(OC_Train),...
    'BasisFunction','constant',...
    'KernelFunction',Hybrid_two,...
    'KernelParameters',theta0,...
    'OptimizeHyperparameters','auto', ...
    'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus',...
    'Kfold',10))
tic
yfit_Train_gpr_Hybrid_two = predict(Mdl_gpr_Hybrid_two,Spectra_Train_FD);
yfit_Test_gpr_Hybrid_two = predict(Mdl_gpr_Hybrid_two,Spectra_Test_FD);
toc
% 性能评估
PI_Train_gpr_Hybrid_two = ModelAssessment...
    (yfit_Train_gpr_Hybrid_two,log(OC_Train),n_Train,1);
PI_Test_gpr_Hybrid_two = ModelAssessment...
    (yfit_Test_gpr_Hybrid_two,log(OC_Test),n_Test,1);

%% 6. 基于不同回归模型
% 支持向量回归--------------------------------------------------
rng default
Mdl_svm = fitrsvm(Spectra_Train_FD,log(OC_Train),...
    'OptimizeHyperparameters','all',...
    'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'))
Mdl_svm.ConvergenceInfo.Converged
% 对训练集和测试集进行预测
tic
yfit_Train_svm = predict(Mdl_svm,Spectra_Train_FD);
yfit_Test_svm = predict(Mdl_svm,Spectra_Test_FD);
toc
% 性能评估
PI_Train_svm = ModelAssessment(yfit_Train_svm,log(OC_Train),n_Train,1);
PI_Test_svm = ModelAssessment(yfit_Test_svm,log(OC_Test),n_Test,1);

% 偏最小二乘回归---------------------------------------------------
LV_max = 30;
LV_PLSR = SelectNumberofLV(Spectra_Train_FD,OC_Train,LV_max);
% LV_num = floor((LV_PLSR.optLV + LV_PLSR.opt_1SD)/2);
tic
[~,~,~,~,beta_pls,PCTvar,mse_pls] = plsregress(Spectra_Train_FD,...
    log(OC_Train),10,'cv',10); % 通过测试不同的验证集结果，LV为19取得最佳结果
toc
tic
yfit_Train_plsr = [ones(n_Train,1) Spectra_Train_FD]*beta_pls;
yfit_Test_plsr = [ones(n_Test,1) Spectra_Test_FD]*beta_pls;
toc
PI_Train_plsr = ModelAssessment(yfit_Train_plsr,log(OC_Train),n_Train,1);
PI_Test_plsr = ModelAssessment(yfit_Test_plsr,log(OC_Test),n_Test,1);

% 随机森林-----------------------------------------------------------
rng default
Mdl_rf = fitrensemble(Spectra_Train_FD,log(OC_Train),...
    'Method','Bag','Learners','tree',...
    'OptimizeHyperparameters',{'NumLearningCycles','MinLeafSize',...
    'MaxNumSplits','NumVariablesToSample'},...
    'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'));
% 预测
tic
yfit_Train_rf = predict(Mdl_rf, Spectra_Train_FD);
yfit_Test_rf = predict(Mdl_rf, Spectra_Test_FD);
toc
% 模型性能评估
PI_Train_rf = ModelAssessment(yfit_Train_rf,log(OC_Train),n_Train,1);
PI_Test_rf = ModelAssessment(yfit_Test_rf,log(OC_Test),n_Test,1);

%% 7.绘制图形
% 图2. 不同SOC含量的光谱图-------------------------------------------------
printSpectra(Spectra_Train,OC_Train,Spectra_Test,OC_Test,0);

printSpectra(Spectra_Train_FD,OC_Train,Spectra_Test_FD,OC_Test,1);
text(500,0,'(b)');

% 图3.不同核函数的观测值与预测值--------------------------------------------
print_OvsP(OC_Test,exp(yfit_Test_gpr_squared_1),150)
text(0,150,'(a)');

print_OvsP(OC_Test,exp(yfit_Test_gpr_rational),150)
text(0,150,'(b)');

print_OvsP(OC_Test,exp(yfit_Test_gpr_ex),150)
text(0,150,'(c)');

print_OvsP(OC_Test,exp(yfit_Test_gpr_matern32),150)
text(0,150,'(d)');

print_OvsP(OC_Test,exp(yfit_Test_gpr_Hybrid_two),150)
text(0,150,'(e)');

% 图4. 不同回归模型的观测值与预测值------------------------------------------
print_OvsP(OC_Test,exp(yfit_Test_plsr),150)
text(0,150,'(a)');

print_OvsP(OC_Test,exp(yfit_Test_svm),150)
text(0,150,'(b)');

print_OvsP(OC_Test,exp(yfit_Test_rf),150)
text(0,150,'(c)');

print_OvsP(OC_Test,exp(yfit_Test_gpr_Hybrid_two),150)
text(0,150,'(d)');

