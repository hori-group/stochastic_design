clear
%% Solve Mode   On:1 Off:0
Analysis=1;
    SolveAve=1; % Compute the upper and lower bounds of Average.
    SolveVar=0; % Compute the upper and lower bounds of Variance. Note that SolveAve must be also 1 if SolveVar is 1.
    SolveCV=1; % Compute the upper and lower bounds of CV. Note that SolveAve must be also 1 if SolveCV is 1.

%% define parameters
n=4; % The number of molecular species in the "repressor" part of the circuit
n2=3; % The number of molecular species in the "reporter" part of the circuit
mu=8; % truncation order
i1=0; % The number of wi=ki reactions in the "repressor" part
i2=6; % The number of wi=kix reactions
i3=0; % The number of wi=kix(x-1) reactions
i4=2; % The number of wi=kix1x2 reactions
j1=0; % The number of wi=ki reactions in the "reporter" part
j2=4; % The number of wi=kix reactions
j3=0; % The number of wi=kix(x-1) reactions
j4=0; % The number of wi=kix1x2 reactions 
Num=20;
param=[0.2 log(2)/5 0.5 log(2)/20 1 1 5 5]; % paramters (repressor)
param2=[0.2 log(2)/5 0.5 log(2)/20]; % parameters (reporter)
trsl=0.1:0.1:0.1+0.1*(Num-1); % translation rates to analyze
halfT=10:10:10+10*(Num-1); % half-life times to analyze
de=zeros(1,Num);
for i=1:Num
    de(1,i)=log(2)/halfT(1,i);
end
FBgain=[5 1 0.1];

X=[0 1 1 0 0 0 0 0;
   0 0 0 1 0 0 1 1;
   1 0 0 0 1 0 1 0;
   0 0 0 0 0 1 0 1]; % Variables of each reaction. (repressor)
X2=[1 0 0 0;
    0 1 1 0;
    0 0 0 1]; % Variables of each reaction (reporter)
S=[1 -1 0 0 0 0 0 0;
   0 0 1 -1 1 1 -1 -1;
   0 0 0 0 1 0 -1 0;
   0 0 0 0 0 1 0 -1]; % Stoichiometry of each reaction (repressor)
S2=[0 0 0 0;
    1 -1 0 0;
    0 0 1 -1]; % Stoichiometry of each reaction (reporter)
neg=[0 0 0 0 1 1 0 0]; % 1 if variable's coefficient is negative in transition rate wi(x). (repressor)
neg2=[0 0 0 0]; % 1 if variable's coefficient is negative in transition rate wi(x). (reporter)
Const=[0 0 0 0 50 50 0 0]; % constants for transition rates wi(x)'s (repressor)
Const2=[0 0 0 0]; % constants for transition rates wi(x)'s (reporter)
Norm=[10 50 10 10]; % Scaling factor for each variable (repressor)
Norm2=[10 20 50]; % Scaling factor for each variable (reporter)

%% Select molecules to analyze
AnaMol=1;
AnaMol2=3;

%% memory allocation
Ave=zeros(Num,Num,2*length(FBgain));
Var=zeros(Num,Num,2*length(FBgain));
CV=zeros(Num,Num,2*length(FBgain));

K=1;
param(1,7)=FBgain(1,K);
param(1,8)=FBgain(1,K);
for i=1:10
    j=8; % Read j-th row entries from the file.
%   for j=1:Num
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
                Ave(i,j,1+2*(K-1))=avemax
                Ave(i,j,2+2*(K-1))=avemin
            end
            if SolveCV==1
                CVmax=OrdervsCVmaxSDP(A,H,mu,n,AnaMol,Arow,Norm);
                CV(i,j,1+2*(K-1))=CVmax
                CVmin=OrdervsCVminSDP(A,H,mu,n,AnaMol,Arow,Norm,avemax);
                CV(i,j,2+2*(K-1))=CVmin
            end

        end
%   end
end
save  SteadyStateAnalysis_TrslDeg_mRNA1.mat