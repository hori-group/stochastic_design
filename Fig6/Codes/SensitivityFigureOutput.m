clear
set(0, 'defaultTextFontName', 'Calibri');
set(0, 'defaultAxesFontName', 'Calibri');
set(0, 'DefaultTextFontSize',20);
set(0, 'defaultAxesFontSize', 20);

load SensitivityAnalysis.mat

figure(1)
rectangle('Position',[0,15,8.5,5],'FaceColor',[1 0.9 0.9])
hold on
for i=1:4
rectangle('Position',[i-0.3,15,0.6,Ave(1,i,1)-15],'FaceColor',[0 1 1])
hold on
rectangle('Position',[i-0.3,15,0.6,Ave(2,i,1)-15],'FaceColor',[0.9 1 1])
hold on
plot([i-0.3,i+0.3],[ave0max(1),ave0max(1)],'b','LineWidth',3)
hold on
rectangle('Position',[i-0.3,15,0.6,5],'FaceColor',[0.9 0.9 0.9])
hold on
end

i=5;
rectangle('Position',[i-0.3+0.5,15,0.6,Ave(1,i,1)-15],'FaceColor',[0 1 1])
hold on
rectangle('Position',[i-0.3+0.5,15,0.6,Ave(2,i,1)-15],'FaceColor',[0.9 1 1])
hold on
plot([i-0.3+0.5,i+0.3+0.5],[ave0max(1),ave0max(1)],'b','LineWidth',3)
hold on
rectangle('Position',[i-0.3+0.5,15,0.6,5],'FaceColor',[0.9 0.9 0.9])
hold on


plot([0,8],[20,20],'r','LineWidth',2)
xlim([0,8])
ylim([15,30])
ylabel('Mean copy number of repressor')

NominalError=(ave0max(1)-ave0min(1))/ave0min(1)*100
AveError=(Ave(:,:,1)-ones(2,5).*((ave0max(1)+ave0min(1))/2))./((ave0max(1)+ave0min(1))/2).*100


figure(2)
rectangle('Position',[0,0.3,8.5,0.1],'FaceColor',[1 0.9 0.9])
hold on
for i=1:4
rectangle('Position',[i-0.3,0.2,0.6,CV(1,i,1)-0.2],'FaceColor',[0 1 1])
hold on
rectangle('Position',[i-0.3,0.2,0.6,CV0max(1)-0.2],'FaceColor',[0.9 1 1])
hold on
plot([i-0.3,i+0.3],[CV0max(1),CV0max(1)],'b','LineWidth',3)
hold on
end
i=5;
rectangle('Position',[i-0.3+0.5,0.2,0.6,CV(1,i,1)-0.2],'FaceColor',[0 1 1])
hold on
rectangle('Position',[i-0.3+0.5,0.2,0.6,CV0max(1)-0.2],'FaceColor',[0.9 1 1])
hold on
plot([i-0.3+0.5,i+0.3+0.5],[CV0max(1),CV0max(1)],'b','LineWidth',3)
hold on
plot([0,8],[0.30,0.30],'r','LineWidth',2)
xlim([0,8])
ylim([0.20,0.35])
ylabel('CV of repressor')