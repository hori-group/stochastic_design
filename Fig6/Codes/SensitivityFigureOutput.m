clear
set(0, 'defaultTextFontName', 'Calibri');
set(0, 'defaultAxesFontName', 'Calibri');
set(0, 'DefaultTextFontSize',20);
set(0, 'defaultAxesFontSize', 20);

%% Load .mat files

load SteadyStateAnalysis_ParamMapSensitivity_1.mat
AveALL=Ave;
CVALL=CV;

load SteadyStateAnalysis_ParamMapSensitivity_2.mat
AveALL=AveALL+Ave;
CVALL=CVALL+CV;

load SteadyStateAnalysis_ParamMapSensitivity_3.mat
AveALL=AveALL+Ave;
CVALL=CVALL+CV;

load SteadyStateAnalysis_ParamMapSensitivity_4.mat
AveALL=AveALL+Ave;
CVALL=CVALL+CV;

load SteadyStateAnalysis_ParamMapSensitivity_5.mat
AveALL=AveALL+Ave;
CVALL=CVALL+CV;

%% Preparation

AveCenter=(AveALL(1,:,:)+AveALL(2,:,:))/2;
AveError=(AveALL(1,:,:)-AveCenter)./AveCenter*100;
maxAveError=max(vec(AveError))
CVCenter=(CVALL(1,:,:)+CVALL(2,:,:))/2;
CVError=(CVALL(1,:,:)-CVCenter)./CVCenter*100;
maxCVError=max(vec(CVError))

AveMAX=permute(AveALL(1,:,:),[3 2 1]);
AveMIN=permute(AveALL(2,:,:),[3 2 1]);
CVMAX=permute(CVALL(1,:,:),[3 2 1]);
CVMIN=permute(CVALL(2,:,:),[3 2 1]);

AveMAX20=AveMAX;
AveMIN20=AveMIN;
CVMAX20=CVMAX;
CVMIN20=CVMIN;

AveMAX10=AveMAX(:,11:31);
AveMIN10=AveMIN(:,11:31);
CVMAX10=CVMAX(:,11:31);
CVMIN10=CVMIN(:,11:31);

AveMAX5=AveMAX(:,16:26);
AveMIN5=AveMIN(:,16:26);
CVMAX5=CVMAX(:,16:26);
CVMIN5=CVMIN(:,16:26);

AveNominal=(AveMAX(1,21)+AveMIN(1,21))/2;
CVNominal=(CVMAX(1,21)+CVMIN(1,21))/2;

%% 20% perturbation

PerMAXAve=max(AveMAX20,[],2)
PerMINAve=-max(-AveMAX20,[],2)
PerMAXCV=max(CVMAX20,[],2)
PerMINCV=-max(-CVMAX20,[],2)


figure(1)
for i=1:4
rectangle('Position',[i-0.3,15,0.6,PerMAXAve(i)-15],'FaceColor',[0.7 1 1])
hold on
rectangle('Position',[i-0.3,15,0.6,PerMINAve(i)-15],'FaceColor',[0.9 1 1])
hold on
rectangle('Position',[i-0.3,15,0.6,5],'FaceColor',[0.9 0.9 0.9])
hold on
end

i=5;
rectangle('Position',[i-0.3+0.5,15,0.6,PerMAXAve(i)-15],'FaceColor',[0.7 1 1])
hold on
rectangle('Position',[i-0.3+0.5,15,0.6,PerMINAve(i)-15],'FaceColor',[0.9 1 1])
hold on
rectangle('Position',[i-0.3+0.5,15,0.6,5],'FaceColor',[0.9 0.9 0.9])
hold on

figure(2)
rectangle('Position',[0,0.3,8.5,0.1],'FaceColor',[1 0.9 0.9])
hold on
for i=1:4
rectangle('Position',[i-0.3,0.2,0.6,PerMAXCV(i)-0.2],'FaceColor',[0.7 1 1])
hold on
rectangle('Position',[i-0.3,0.2,0.6,PerMINCV(i)-0.2],'FaceColor',[0.9 1 1])
hold on
end
i=5;
rectangle('Position',[i-0.3+0.5,0.2,0.6,PerMAXCV(i)-0.2],'FaceColor',[0.7 1 1])
hold on
rectangle('Position',[i-0.3+0.5,0.2,0.6,PerMINCV(i)-0.2],'FaceColor',[0.9 1 1])
hold on


%% 10% perturbation

PerMAXAve=max(AveMAX10,[],2)
PerMINAve=-max(-AveMAX10,[],2)
PerMAXCV=max(CVMAX10,[],2)
PerMINCV=-max(-CVMAX10,[],2)

figure(1)
rectangle('Position',[0,15,8.5,5],'FaceColor',[1 0.9 0.9])
hold on
for i=1:4
rectangle('Position',[i-0.2,PerMINAve(i),0.4,PerMAXAve(i)-PerMINAve(i)],'FaceColor',[0.4 0.93 0.93])
hold on
end

i=5;
rectangle('Position',[i-0.2+0.5,PerMINAve(i),0.4,PerMAXAve(i)-PerMINAve(i)],'FaceColor',[0.4 0.93 0.93])
hold on

figure(2)
for i=1:4
rectangle('Position',[i-0.2,PerMINCV(i),0.4,PerMAXCV(i)-PerMINCV(i)],'FaceColor',[0.4 0.93 0.93])
hold on
end
i=5;
rectangle('Position',[i-0.2+0.5,PerMINCV(i),0.4,PerMAXCV(i)-PerMINCV(i)],'FaceColor',[0.4 0.93 0.93])
hold on

%% 5% perturbation

PerMAXAve=max(AveMAX5,[],2)
PerMINAve=-max(-AveMAX5,[],2)
PerMAXCV=max(CVMAX5,[],2)
PerMINCV=-max(-CVMAX5,[],2)

figure(1)
rectangle('Position',[0,15,8.5,5],'FaceColor',[1 0.9 0.9])
hold on
for i=1:4
rectangle('Position',[i-0.1,PerMINAve(i),0.2,PerMAXAve(i)-PerMINAve(i)],'FaceColor',[0 0.85 0.85])
hold on
plot([i-0.3,i+0.3],[AveNominal,AveNominal],'b','LineWidth',3)
hold on
rectangle('Position',[i-0.3,15,0.6,5],'FaceColor',[0.9 0.9 0.9])
hold on
end

i=5;
rectangle('Position',[i-0.1+0.5,PerMINAve(i),0.2,PerMAXAve(i)-PerMINAve(i)],'FaceColor',[0 0.85 0.85])
hold on
plot([i-0.3+0.5,i+0.3+0.5],[AveNominal,AveNominal],'b','LineWidth',3)
hold on
rectangle('Position',[i-0.3+0.5,15,0.6,5],'FaceColor',[0.9 0.9 0.9])
hold on


plot([0,8],[20,20],'r','LineWidth',2)
xlim([0,8])
ylim([15,30])
ylabel('Mean copy number of repressor')


figure(2)
for i=1:4
rectangle('Position',[i-0.1,PerMINCV(i),0.2,PerMAXCV(i)-PerMINCV(i)],'FaceColor',[0 0.85 0.85])
hold on
plot([i-0.3,i+0.3],[CVNominal,CVNominal],'b','LineWidth',3)
hold on
end
i=5;
rectangle('Position',[i-0.1+0.5,PerMINCV(i),0.2,PerMAXCV(i)-PerMINCV(i)],'FaceColor',[0 0.85 0.85])
hold on
plot([i-0.3+0.5,i+0.3+0.5],[CVNominal,CVNominal],'b','LineWidth',3)
hold on
plot([0,8],[0.30,0.30],'r','LineWidth',2)
xlim([0,8])
ylim([0.20,0.35])
ylabel('CV of repressor')