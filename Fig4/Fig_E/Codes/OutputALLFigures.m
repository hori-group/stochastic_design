clear;
load('SteadyStateAnalysis_TrslDeg_Reporter_CV_1.mat');
AveALL=Ave;
CVALL=CV;
load('SteadyStateAnalysis_TrslDeg_Reporter_CV_2.mat');
CVALL=CV+CVALL;
load('SteadyStateAnalysis_TrslDeg_Reporter_CV_3.mat');
CVALL=CV+CVALL;
load('SteadyStateAnalysis_TrslDeg_Reporter_CV_4.mat');
CVALL=CV+CVALL;
load('SteadyStateAnalysis_TrslDeg_Reporter_CV_5_1.mat');
CVALL=CV+CVALL;
load('SteadyStateAnalysis_TrslDeg_Reporter_CV_5_2.mat');
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