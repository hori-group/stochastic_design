AveALL=zeros(40,6);
CVALL=zeros(40,6);

load SteadyStateAnalysis_AvevsCV_Repressor_K1.mat
AveALL(:,1)=Ave(:,:,1);
AveALL(:,2)=Ave(:,:,2);
CVALL(:,1)=CV(:,:,1);
CVALL(:,2)=CV(:,:,2);

load SteadyStateAnalysis_AvevsCV_Repressor_K2.mat
AveALL(:,3)=Ave(:,:,3);
AveALL(:,4)=Ave(:,:,4);
CVALL(:,3)=CV(:,:,3);
CVALL(:,4)=CV(:,:,4);

load SteadyStateAnalysis_AvevsCV_Repressor_K3.mat
AveALL(:,5)=Ave(:,:,5);
AveALL(:,6)=Ave(:,:,6);
CVALL(:,5)=CV(:,:,5);
CVALL(:,6)=CV(:,:,6);

for K=1:length(FBgain)
    
    figure(2*length(FBgain)+1)
    plot((AveALL(:,1+2*(K-1))+AveALL(:,2+2*(K-1)))/2,CVALL(:,1+2*(K-1)),'-')
    hold on
    plot((AveALL(:,1+2*(K-1))+AveALL(:,2+2*(K-1)))/2,CVALL(:,2+2*(K-1)),'-')
    hold on
    xlabel('Average')
    ylabel('Coefficient of variation')
end
xlim([0 200])