clear
load SteadyStateAnalysis_RepressorReporter_Relation.mat 
AveTetRALL=AveTetR;
AveGFPALL=AveGFP;

%% Compute copy number with deterministic model
x2=((AveTetRALL(:,1)+AveTetRALL(:,2))./2)';
DT=50;
k7=0.2;
k8=log(2)/5;
k9=0.5
k10=log(2)/20
x6=DT*k7*k9/k8/k10./(1+5*x2);

%% Output

%yyaxis left
plot((AveTetRALL(:,1)+AveTetRALL(:,2))/2,AveGFPALL(:,1),'--o');
hold on
%yyaxis right
plot(x2,x6,'--x');
legend('Stochastic model','Deterministic model')

AveGFPError=(AveGFPALL(:,1)-AveGFPALL(:,2))./AveGFPALL(:,1)*100 % Evaluation of difference between the upper and lower bounds (%). 
DetSDPError2=(x6'-AveGFPALL(:,2))./AveGFPALL(:,1)*100 % Evaluation of difference between the result of our rigorous approach (SDP) and that of deterministic model (%). 