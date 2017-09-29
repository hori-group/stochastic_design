clear;
load('SteadyStateAnalysis_TrslDeg_Repressor_1.mat');
AveALL=Ave;
CVALL=CV;
load('SteadyStateAnalysis_TrslDeg_Repressor_2.mat');
AveALL=Ave+AveALL;
CVALL=CV+CVALL;
load('SteadyStateAnalysis_TrslDeg_Repressor_3.mat');
AveALL=Ave+AveALL;
CVALL=CV+CVALL;
load('SteadyStateAnalysis_TrslDeg_Repressor_4.mat');
AveALL=Ave+AveALL;
CVALL=CV+CVALL;
load('SteadyStateAnalysis_TrslDeg_Repressor_5.mat');
AveALL=Ave+AveALL;
CVALL=CV+CVALL;

figure(1)
pcolor(halfT,trsl,AveALL(:,:,2));
xlabel('Half degradation time')
ylabel('Translation gain')
c=colorbar;
c.Label.String = 'Mean copy number'
figure(2)
pcolor(halfT,trsl,CVALL(:,:,1));
xlabel('Half degradation time')
ylabel('Translation gain')
c=colorbar;
c.Label.String = 'Coefficient of variation'