% 绘制不同范围的
function printSpectra(Spectra_Train,OC_Train,Spectra_Test,OC_Test,order)

if order == 0
    wavelengths = (500:2:2498);
else
    wavelengths = (500:2:2496);
end

% 不同含量范围的有机碳的索引号
idx_0to20_Train = find(0<=OC_Train & OC_Train<=20);
idx_20to40_Train = find(20<OC_Train & OC_Train<=40);
idx_40to60_Train = find(40<OC_Train & OC_Train<=60);
idx_60to75_Train = find(60<OC_Train & OC_Train<=75);
idx_75to125_Train = find(75<OC_Train & OC_Train<=125);
idx_125to200_Train = find(125<OC_Train & OC_Train<=200);
idx_200to350_Train = find(200<OC_Train & OC_Train<=350);
idx_350to550_Train = find(350<OC_Train & OC_Train<=550);

% 不同含量范围的有机碳的索引号
idx_0to20_Test = find(0<=OC_Test & OC_Test<=20);
idx_20to40_Test = find(20<OC_Test & OC_Test<=40);
idx_40to60_Test = find(40<OC_Test & OC_Test<=60);
idx_60to75_Test = find(60<OC_Test & OC_Test<=75);
idx_75to125_Test = find(75<OC_Test & OC_Test<=125);
idx_125to200_Test = find(125<OC_Test & OC_Test<=200);
idx_200to350_Test = find(200<OC_Test & OC_Test<=350);
idx_350to550_Test = find(350<OC_Test & OC_Test<=550);

figure;
plot(wavelengths,...
    [mean([Spectra_Train(idx_0to20_Train,:);Spectra_Test(idx_0to20_Test,:)]);...
    mean([Spectra_Train(idx_20to40_Train,:);Spectra_Test(idx_20to40_Test,:)]);...
    mean([Spectra_Train(idx_40to60_Train,:);Spectra_Test(idx_40to60_Test,:)]);...
    mean([Spectra_Train(idx_60to75_Train,:);Spectra_Test(idx_60to75_Test,:)]);...
    mean([Spectra_Train(idx_75to125_Train,:);Spectra_Test(idx_75to125_Test,:)]);...
    mean([Spectra_Train(idx_125to200_Train,:);Spectra_Test(idx_125to200_Test,:)]);...
    mean([Spectra_Train(idx_200to350_Train,:);Spectra_Test(idx_200to350_Test,:)]);...
    mean([Spectra_Train(idx_350to550_Train,:);Spectra_Test(idx_350to550_Test,:)])],...
    'Linewidth',1);
xlim([400 2600]); % 设置 x轴的范围
xlabel('Wavelength (nm)');
ylabel('Absorbance');
grid on
grid minor % 小网格线
lg1 = legend('1-20','20-40','40-60','60-75','75-125',...
    '125-200','200-350','350-550');
title(lg1,'SOC Class (g/kg)');
end
