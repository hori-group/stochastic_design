clear;

%% Repressor Load
load('SteadyStateAnalysis_TrslDeg_Repressor_1.mat');
RepressorAve=Ave;
RepressorCV=CV;
load('SteadyStateAnalysis_TrslDeg_Repressor_2.mat');
RepressorAve=Ave+RepressorAve;
RepressorCV=CV+RepressorCV;
load('SteadyStateAnalysis_TrslDeg_Repressor_3.mat');
RepressorAve=Ave+RepressorAve;
RepressorCV=CV+RepressorCV;
load('SteadyStateAnalysis_TrslDeg_Repressor_4.mat');
RepressorAve=Ave+RepressorAve;
RepressorCV=CV+RepressorCV;
load('SteadyStateAnalysis_TrslDeg_Repressor_5.mat');
RepressorAve=Ave+RepressorAve;
RepressorCV=CV+RepressorCV;

%% Reporter Load
load('SteadyStateAnalysis_TrslDeg_Reporter_CV_1.mat');
ReporterAve=Ave;
ReporterCV=CV;
load('SteadyStateAnalysis_TrslDeg_Reporter_CV_2.mat');
ReporterAve=Ave+ReporterAve;
ReporterCV=CV+ReporterCV;
load('SteadyStateAnalysis_TrslDeg_Reporter_CV_3.mat');
ReporterAve=Ave+ReporterAve;
ReporterCV=CV+ReporterCV;
load('SteadyStateAnalysis_TrslDeg_Reporter_CV_4.mat');
ReporterAve=Ave+ReporterAve;
ReporterCV=CV+ReporterCV;
load('SteadyStateAnalysis_TrslDeg_Reporter_CV_5_1.mat');
ReporterAve=Ave+ReporterAve;
ReporterCV=CV+ReporterCV;
load('SteadyStateAnalysis_TrslDeg_Reporter_CV_5_2.mat');
ReporterAve=Ave+ReporterAve;
ReporterCV=CV+ReporterCV;


X=zeros(size(RepressorCV,1),size(RepressorCV,2));

for i=1:size(RepressorCV,1)
    for j=1:size(RepressorCV,2)
        if ReporterCV(i,j,1)<=0.9 && RepressorCV(i,j,1)<=0.30 && RepressorAve(i,j,2)>=20
            X(i,j)=1;
        end
    end
end

figure(1)
pcolor(halfT,trsl,X);
xlabel('Half degradation time')
ylabel('Translation gain')
