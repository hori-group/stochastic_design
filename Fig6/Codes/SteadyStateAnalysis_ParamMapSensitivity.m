clear
%% Solve Mode   On:1 Off:0
Analysis=1;
    SolveAve=1; % This parameter computes the upper and lower bounds of Average.
    SolveVar=0; % This parameter computes the upper and lower bounds of Variance. Note that SolveAve must be also 1 if SolveVar is 1.
    SolveCV=1; % This parameter computes the upper and lower bounds of CV. Note that SolveAve must be also 1 if SolveCV is 1.

%% define parameters
n=4; % The number of variables (molecular species)
mu=8; % The truncation order of the system
i1=0; % The number of wi=ki reactions
i2=6; % The number of wi=kix reactions
i3=0; % The number of wi=kix(x-1) reactions
i4=2; % The number of wi=kix1x2 reactions
param0=[0.2 log(2)/5 0.4 log(2)/80 1 1 5 5]; % Parameters of reactions. Refer "Readme.txt" for details.
PerRate=-0.2:0.01:0.2; % The rate of perturbation
Num=length(PerRate);

X=[0 1 1 0 0 0 0 0;
   0 0 0 1 0 0 1 1;
   1 0 0 0 1 0 1 0;
   0 0 0 0 0 1 0 1]; % Variables of each reaction. Refer "Readme.txt" for details. 
S=[1 -1 0 0 0 0 0 0;
   0 0 1 -1 1 1 -1 -1;
   0 0 0 0 1 0 -1 0;
   0 0 0 0 0 1 0 -1]; % Stoichiometry of each reaction. Refer "Readme.txt" for details.
neg=[0 0 0 0 1 1 0 0]; % Put 1 if variable's coefficient is negative in transition rate wi(x). Refer "Readme.txt" for details.
Const0=[0 0 0 0 50 50 0 0]; % Put the number of constant if transition rate wi(x)'s variable contains constant terms. Refer "Readme.txt" for details.
Norm=[20 50 10 10]; % Scaling number for each variable

%% Choose Molecule to analyze
AnaMol=2;

%% Preparation (Memory Allocation)
ParamMap=zeros(length(param0),Num,5);
ConstMap=zeros(length(Const0),Num,5);
Ave=zeros(2,Num,5);
Var=zeros(2,Num,5);
CV=zeros(2,Num,5);

%% Calculation

for i=1:5
    for j=1:Num
        %% Parameter
        
        param=param0;
        Const=Const0;
        
        if i<=4
            param(i)=param0(i)*(1+PerRate(j));
        else
            Const(5)=Const0(5)*(1+PerRate(j));
        end
        ParamMap(:,j,i)=param';
        ConstMap(:,j,i)=Const';

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
                Ave(1,j,i)=avemax
                Ave(2,j,i)=avemin
            end
            if SolveVar==1
                Varmax=OrdervsVarmaxSDP(A,H,mu,n,AnaMol,Arow,Norm);
                Var(1,j,i)=Varmax
                Varmin=OrdervsVarminSDP(A,H,mu,n,AnaMol,Arow,Norm,avemax);
                Var(2,j,i)=Varmin
            end
            if SolveCV==1
                CVmax=OrdervsCVmaxSDP(A,H,mu,n,AnaMol,Arow,Norm);
                CV(1,j,i)=CVmax
                CVmin=OrdervsCVminSDP(A,H,mu,n,AnaMol,Arow,Norm,avemax);
                CV(2,j,i)=CVmin
            end

        end
    end
end

save SteadyStateAnalysis_ParamMapSensitivity.mat