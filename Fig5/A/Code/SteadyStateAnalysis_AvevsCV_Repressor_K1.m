clear
%% Solve Mode   On:1 Off:0
Analysis=1;
    SolveAve=1; % Compute the upper and lower bounds of Average.
    SolveCV=1; % Compute the upper and lower bounds of CV. Note that SolveAve must be also 1 if SolveCV is 1.

%% define parameters
n=4; % The number of variables (molecular species)
mu=8; % Truncation order
i1=0; % The number of wi=ki reactions
i2=6; % The number of wi=kix reactions
i3=0; % The number of wi=kix(x-1) reactions
i4=2; % The number of wi=kix1x2 reactions
Num=40; % The number of plots for analysis
param=[0.2 log(2)/5 0.4 log(2)/80 1 1 10 10]; % reaction parameters
trsl=0.025:0.025:0.025+0.025*(Num-1);
FBgain=[5 1 0.1];

X=[0 1 1 0 0 0 0 0;
   0 0 0 1 0 0 1 1;
   1 0 0 0 1 0 1 0;
   0 0 0 0 0 1 0 1]; % Variables of each reaction. 
S=[1 -1 0 0 0 0 0 0;
   0 0 1 -1 1 1 -1 -1;
   0 0 0 0 1 0 -1 0;
   0 0 0 0 0 1 0 -1]; % Stoichiometry of each reaction.
neg=[0 0 0 0 1 1 0 0]; % 1 if variable's coefficient is negative in transition rate wi(x). 
Const=[0 0 0 0 50 50 0 0]; % Constants for transition rates wi(x)'s
Norm=[10 20 10 10]; % Scaling factor for each variable

%% Select molecule to analyze
AnaMol=2;

%% memory allocation
Ave=zeros(Num,1,2*length(FBgain));
Var=zeros(Num,1,2*length(FBgain));
CV=zeros(Num,1,2*length(FBgain));

K=1;
param(1,7)=FBgain(1,K);
param(1,8)=FBgain(1,K);
for i=1:Num
    %% Parameter
    param(1,3)=trsl(1,i);

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
            Ave(i,1,1+2*(K-1))=avemax
            Ave(i,1,2+2*(K-1))=avemin
        end
        if SolveCV==1
            CVmax=OrdervsCVmaxSDP(A,H,mu,n,AnaMol,Arow,Norm);
            CV(i,1,1+2*(K-1))=CVmax
            CVmin=OrdervsCVminSDP(A,H,mu,n,AnaMol,Arow,Norm,avemax);
            CV(i,1,2+2*(K-1))=CVmin
        end

    end
end
save SteadyStateAnalysis_AvevsCV_Repressor_K1.mat