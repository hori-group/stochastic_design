clear
%% Solve Mode   On:1 Off:0
Analysis=1;
    SolveAve=1; % Compute the upper and lower bounds of Average.
    SolveVar=0; % Compute the upper and lower bounds of Variance. Note that SolveAve must be also 1 if SolveVar is 1.
    SolveCV=1; % Compute the upper and lower bounds of CV. Note that SolveAve must be also 1 if SolveCV is 1.

%% define parameters
n=4; % The number of variables (molecular species)
mu=8; % The truncation order of the system
i1=0; % The number of wi=ki reactions (zero-th order)
i2=6; % The number of wi=kix reactions (first order)
i3=0; % The number of wi=kix(x-1) reactions (second order)
i4=2; % The number of wi=kix1x2 reactions (second order)
Num=20; % The number of plots for each parameter 
param=[0.2 log(2)/5 0.5 log(2)/20 1 1 5 5]; % reaction parameters
trsl=0.05:0.05:0.05+0.05*(Num-1); % set of translation rates to analyze
halfT=5:5:5+5*(Num-1); % set of the half-life times of protein to analyze

de=zeros(1,Num);
for i=1:Num
    de(1,i)=log(2)/halfT(1,i); % Transform half-life times into degradation rates
end

X=[0 1 1 0 0 0 0 0;
   0 0 0 1 0 0 1 1;
   1 0 0 0 1 0 1 0;
   0 0 0 0 0 1 0 1]; % Variables of each reaction. 
S=[1 -1 0 0 0 0 0 0;
   0 0 1 -1 1 1 -1 -1;
   0 0 0 0 1 0 -1 0;
   0 0 0 0 0 1 0 -1]; % Stoichiometry of each reaction. 
neg=[0 0 0 0 1 1 0 0]; % 1 if variable's coefficient is negative in transition rate wi(x).
Const=[0 0 0 0 50 50 0 0]; % Constant term of transition rate wi(x). 
Norm=[20 50 10 10]; % Scaling factor for each variable

%% Choose Molecule to analyze
AnaMol=2;

%% Preparation (Memory Allocation)
Ave=zeros(Num,Num,2);
Var=zeros(Num,Num,2);
CV=zeros(Num,Num,2);

%% Calculation

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
        A=Amatrix(mu,n,I,param,X,S,neg,Const,Acolumn,Arow);
        A=ANorm(A,Norm,Acolumn,Arow,mu,n); % Normalization

        %% Stieltjes Moment Condition

        H=MomentCondition(mu,n,Arow);

        %% Solve SDP
        if Analysis==1
            if SolveAve==1
                [avemax,avemin]=OrdervsAveSDP(A,H,mu,n,AnaMol,Arow,Norm);
                Ave(i,j,1)=avemax
                Ave(i,j,2)=avemin
            end
            if SolveVar==1
                Varmax=OrdervsVarmaxSDP(A,H,mu,n,AnaMol,Arow,Norm);
                Var(i,j,1)=Varmax
                Varmin=OrdervsVarminSDP(A,H,mu,n,AnaMol,Arow,Norm,avemax);
                Var(i,j,2)=Varmin
            end
            if SolveCV==1
                CVmax=OrdervsCVmaxSDP(A,H,mu,n,AnaMol,Arow,Norm);
                CV(i,j,1)=CVmax
                CVmin=OrdervsCVminSDP(A,H,mu,n,AnaMol,Arow,Norm,avemax);
                CV(i,j,2)=CVmin
            end

        end
    end
end

save SteadyStateAnalysis_TrslDeg_Repressor.mat

