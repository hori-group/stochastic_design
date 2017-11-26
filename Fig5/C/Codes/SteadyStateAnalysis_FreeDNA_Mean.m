clear
%% Solve Mode   On:1 Off:0
Analysis=1;
    SolveAve=1; % This parameter computes the upper and lower bounds of Average.

%% define parameters
n=4; % The number of molecular species at the reaction of repressing process
mu=8; % The truncation order of the moments.
i1=0; % The number of wi=ki reactions at the reaction of repressing process
i2=6; % The number of wi=kix reactions at the reaction of repressing process
i3=0; % The number of wi=kix(x-1) reactions at the reaction of repressing process
i4=2; % The number of wi=kix1x2 reactions at the reaction of repressing process
param=[0.2 log(2)/5 0.5 log(2)/80 1 1 5 5]; % Parameters of reaction at the repressing process. Refer "Readme.txt" for details.
trsl=[0.01:0.01:0.04 0.05 0.07 0.10 0.15 0.25 0.50 1.00]; % The vector of all the translation rates to analyze
Num=length(trsl);

X=[0 1 1 0 0 0 0 0;
   0 0 0 1 0 0 1 1;
   1 0 0 0 1 0 1 0;
   0 0 0 0 0 1 0 1]; % Variables of each reaction at the repressing process. Refer "Readme.txt" for details. 
S=[1 -1 0 0 0 0 0 0;
   0 0 1 -1 1 1 -1 -1;
   0 0 0 0 1 0 -1 0;
   0 0 0 0 0 1 0 -1]; % Stoichiometry of each reaction at the repressing process. Refer "Readme.txt" for details.
neg=[0 0 0 0 1 1 0 0]; % Put 1 if variable's coefficient is negative in transition rate wi(x). Refer "Readme.txt" for details. (repressing process)
Const=[0 0 0 0 50 50 0 0]; % Put the number of constant if transition rate wi(x)'s variable contains constant terms. Refer "Readme.txt" for details. (repressing process)
Norm=[10 20 10 10]; % Scaling number for each variable (repressing process)

%% Choose Molecule to analyze
AnaMol=4; % Repressing process

%% Preparation (Memory Allocation)
Ave=zeros(Num,2);

%% Calculation

for i=1:Num
    %% Parameter
    param(1,3)=trsl(i);

    %% Location of moments

    Acolumn=momentlocation(mu,n);
    Arow=[zeros(n,1),momentlocation(mu+1,n)];

    %% Determine A matrix

    I=[i1 i2 i3 i4];
    A=Amatrix(mu,n,I,param,X,S,neg,Const,Acolumn,Arow);
    A=ANorm(A,Norm,Acolumn,Arow,mu,n); % Normalization

    %% Stieltjes Moment Condition

    H=MomentCondition(mu,n,Arow);

    %% Solve SDP
    if Analysis==1
        if SolveAve==1
            [avemax,avemin]=OrdervsAveSDP(A,H,mu,n,AnaMol,Arow,Norm);
            Ave(i,1)=avemax
            Ave(i,2)=avemin
        end    
    end
end
    
save SteadyStateAnalysis_FreeDNA_Mean.mat