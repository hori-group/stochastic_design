clear
load SteadyStateAnalysis_TrslDeg_Reporter_Ave.mat
clearvars -except Ave % memory saving

%% define parameters
n=6; % The number of variables (molecular species)
mu=6; % Truncation order
i1=0; % The number of wi=ki reactions
i2=10; % The number of wi=kix reactions
i3=0; % The number of wi=kix(x-1) reactions
i4=2; % The number of wi=kix1x2 reactions
Num=20;
param=[0.2 log(2)/5 0.5 log(2)/20 1 1 0.2 log(2)/5 0.5 log(2)/20 5 5]; % reaction parameters
trsl=0.05:0.05:0.05+0.05*(Num-1); % Translation rates to analyze
halfT=5:5:5+5*(Num-1); % Half-life times of protein to analyze
de=zeros(1,Num);
for i=1:Num
    de(1,i)=log(2)/halfT(1,i); % Transform half-life times to degradation rates
end

X=[0 1 1 0 0 0 0 0 0 0 0 0;
   0 0 0 1 0 0 0 0 0 0 1 1;
   1 0 0 0 1 0 0 0 0 0 1 0;
   0 0 0 0 0 1 1 0 0 0 0 1;
   0 0 0 0 0 0 0 1 1 0 0 0;
   0 0 0 0 0 0 0 0 0 1 0 0]; % Variables of each reaction. 
S=[1 -1 0 0 0 0 0 0 0 0 0 0;
   0 0 1 -1 1 1 0 0 0 0 -1 -1;
   0 0 0 0 1 0 0 0 0 0 -1 0;
   0 0 0 0 0 1 0 0 0 0 0 -1;
   0 0 0 0 0 0 1 -1 0 0 0 0;
   0 0 0 0 0 0 0 0 1 -1 0 0]; % Stoichiometry of each reaction.
neg=[0 0 0 0 1 1 0 0 0 0 0 0]; % 1 if variable's coefficient is negative in transition rate wi(x). 
Const=[0 0 0 0 50 50 0 0 0 0 0 0]; % constants for transition rate wi(x)'s.
Norm=[10 20 10 10 10 20]; % Scaling factor for each variable

%% Select molecules to analyze
AnaMol=6;

%% memory allocation
CV=zeros(Num,Num,2);

for i=1:Num
    for j=1:Num
        %% Parameter
        param(1,3)=trsl(1,i);
        param(1,4)=de(1,j);

        %% Location of moments

        Acolumn=momentlocation(mu,n);
        Arow=[zeros(n,1),momentlocation(mu+1,n)];

        %% Determine A matrix

        I=[i1 i2 i3 i4];
        A=Amatrix_Fast(mu,n,I,param,X,S,neg,Const,Acolumn,Arow);
        A=ANorm(A,Norm,Acolumn,Arow,mu,n); % Normalization

        %% Stieltjes Moment Condition

        H=MomentCondition_Fast(mu,n,Arow);

        %% Solve SDP
        avemax=Ave(i,j,1);
        avemin=Ave(i,j,2);
        CVmax=OrdervsCVmaxSDP_Reporter(A,H,mu,n,AnaMol,Arow,Norm);
        CV(i,j,1)=CVmax
        CVmin=OrdervsCVminSDP_Reporter(A,H,mu,n,AnaMol,Arow,Norm,avemax);
        CV(i,j,2)=CVmin
     end
 end
         
save SteadyStateAnalysis_TrslDeg_Reporter_CV.mat