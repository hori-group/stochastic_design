clear;

%% Load Files

load('SteadyStateAnalysis_Reporter_Var_k9_05.mat');
ReporVar=Var;

load('SteadyStateAnalysis_Reporter_Var_k9_1.mat');
ReporVar=cat(3,ReporVar,Var);

load('SteadyStateAnalysis_Reporter_Var_k9_2.mat');
ReporVar=cat(3,ReporVar,Var);

load('SteadyStateAnalysis_FreeDNA_Variance.mat')
FreeDNAVar=Var;

ReporVar(7,:,:)=[];
FreeDNAVar(7,:,:)=[];
TranslRate=trsl;
TranslRate(7)=[];

%% Load Files (Add files)

load('SteadyStateAnalysis_Reporter_Var_k9_05_Add.mat');
ReporVarAdd=Var;

load('SteadyStateAnalysis_Reporter_Var_k9_1_Add.mat');
ReporVarAdd=cat(3,ReporVarAdd,Var);

load('SteadyStateAnalysis_Reporter_Var_k9_2_Add.mat');
ReporVarAdd=cat(3,ReporVarAdd,Var);

load('SteadyStateAnalysis_FreeDNA_Variance_Add.mat')
FreeDNAVar=cat(1,Var,FreeDNAVar);

ReporVar=cat(1,ReporVarAdd,ReporVar);
TranslRate=[trsl,TranslRate];

%% Preparation

ReporVarCenter=permute((ReporVar(:,1,:)+ReporVar(:,2,:))/2,[1 3 2]);
ReporVarDist=permute((ReporVar(:,1,:)-ReporVar(:,2,:))/2,[1 3 2]);
ReporSDPError=ReporVarDist./ReporVarCenter*100
FreeDNAVarCenter=(FreeDNAVar(:,1)+FreeDNAVar(:,2))/2;
FreeDNAVarDist=(FreeDNAVar(:,1)-FreeDNAVar(:,2))/2;
FreeDNASDPError=FreeDNAVarDist./FreeDNAVarCenter*100

%% Output Figures

figure(1)
plot(FreeDNAVarCenter,ReporVarCenter(:,1),'--o');
xlabel('Variance of free DNA')
ylabel('Variance of reporter')
hold on
plot(FreeDNAVarCenter,ReporVarCenter(:,2),'--x');
xlabel('Variance of free DNA')
ylabel('Variance of reporter')
hold on
plot(FreeDNAVarCenter,ReporVarCenter(:,3),'--x');
xlabel('Variance of free DNA')
ylabel('Variance of reporter')
legend('k9=0.5','k9=1.0','k9=2.0')

figure(2)
semilogy(FreeDNAVarCenter,ReporVarCenter(:,1),'--o');
xlabel('Variance of free DNA')
ylabel('Variance of reporter')
hold on
semilogy(FreeDNAVarCenter,ReporVarCenter(:,2),'--x');
xlabel('Variance of free DNA')
ylabel('Variance of reporter')
hold on
semilogy(FreeDNAVarCenter,ReporVarCenter(:,3),'--*');
xlabel('Variance of free DNA')
ylabel('Variance of reporter')
legend('k9=0.5','k9=1.0','k9=2.0')

save Variance_Report_Relation.mat