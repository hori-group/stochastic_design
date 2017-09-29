clear

load SteadyStateAnalysis_TrslDeg_mRNA1.mat

AveALL=zeros(10,2);
CVALL=zeros(10,2);

j=8;
for i=1:10
    AveALL(i,1)=Ave(i,j,1);
    AveALL(i,2)=Ave(i,j,2);
    CVALL(i,1)=CV(i,j,1);
    CVALL(i,2)=CV(i,j,2);
end

trsl=0.1:0.1:1;

yyaxis left
plot(trsl',(AveALL(:,1)+AveALL(:,2))/2,'--x');
yyaxis right
plot(trsl',(CVALL(:,1)+CVALL(:,2))/2,'--o')