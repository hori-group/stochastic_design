clear
load SteadyStateAnalysis_FreeDNA_Mean.mat
clearvars -except Ave trsl Num % Memory Reduction (No effects on results)

%% define parameters
n=6; % The number of variables (molecular species)
mu=6; % The truncation order of the system
i1=0; % The number of wi=ki reactions
i2=10; % The number of wi=kix reactions
i3=0; % The number of wi=kix(x-1) reactions
i4=2; % The number of wi=kix1x2 reactions
param=[0.2 log(2)/5 0.5 log(2)/80 1 1 0.2 log(2)/5 2.0 log(2)/20 5 5]; % Parameters of reaction at the repressing process. Refer "Readme.txt" for details.
param2=[0.2 log(2)/5 2.0 log(2)/20]; % Parameters of reaction at the report process. Refer "Readme.txt" for details.
trsl=[0.01:0.01:0.04 0.05 0.07 0.10 0.15 0.25 0.50 1.00]; % The vector of all the translation rates to analyze

X=[0 1 1 0 0 0 0 0 0 0 0 0;
   0 0 0 1 0 0 0 0 0 0 1 1;
   1 0 0 0 1 0 0 0 0 0 1 0;
   0 0 0 0 0 1 1 0 0 0 0 1;
   0 0 0 0 0 0 0 1 1 0 0 0;
   0 0 0 0 0 0 0 0 0 1 0 0]; % Variables of each reaction. Refer "Readme.txt" for details. 
S=[1 -1 0 0 0 0 0 0 0 0 0 0;
   0 0 1 -1 1 1 0 0 0 0 -1 -1;
   0 0 0 0 1 0 0 0 0 0 -1 0;
   0 0 0 0 0 1 0 0 0 0 0 -1;
   0 0 0 0 0 0 1 -1 0 0 0 0;
   0 0 0 0 0 0 0 0 1 -1 0 0]; % Stoichiometry of each reaction. Refer "Readme.txt" for details.
neg=[0 0 0 0 1 1 0 0 0 0 0 0]; % Put 1 if variable's coefficient is negative in transition rate wi(x). Refer "Readme.txt" for details.
Const=[0 0 0 0 50 50 0 0 0 0 0 0]; % Put the number of constant if transition rate wi(x)'s variable contains constant terms. Refer "Readme.txt" for details.
Norm=[10 20 10 10 10 80]; % Scaling number for each variable

%% Choose Molecule to analyze
AnaMol=6;

%% Preparation (Memory Allocation)
Var=zeros(Num,2);

for i=1:Num
    %% Parameter
    param(1,3)=trsl(i);

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
    avemax=Ave(i,1)*(param2(1,1)*param2(1,3))/(param2(1,2)*param2(1,4));
    avemin=Ave(i,2)*(param2(1,1)*param2(1,3))/(param2(1,2)*param2(1,4));
    Varmax=OrdervsVarmaxSDP(A,H,mu,n,AnaMol,Arow,Norm);
    Var(i,1)=Varmax
    Varmin=OrdervsVarminSDP(A,H,mu,n,AnaMol,Arow,Norm,avemax);
    Var(i,2)=Varmin
 end
         
save SteadyStateAnalysis_Reporter_Var_k9_2.mat