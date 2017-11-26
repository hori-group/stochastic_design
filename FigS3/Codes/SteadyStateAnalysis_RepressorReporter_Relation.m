clear
%% Solve Mode   On:1 Off:0
Analysis=1;
    SolveAve=1; % Compute the upper and lower bounds of Average.
    SolveVar=0; % always 0
    SolveCV=0; % always 0

%% define parameters
n=4; % The number of molecular species in the "repressor" part of the circuit
n2=3; % The number of molecular species in the "reporter" part of the circuit
mu=8; % truncation order of the moments.
i1=0; % The number of wi=ki reactions (repressor)
i2=6; % The number of wi=kix reactions 
i3=0; % The number of wi=kix(x-1) reactions
i4=2; % The number of wi=kix1x2 reactions 
j1=0; % The number of wi=ki reactions (reporter)
j2=4; % The number of wi=kix reactions 
j3=0; % The number of wi=kix(x-1) reactions 
j4=0; % The number of wi=kix1x2 reactions
Num=20; % The number of plots for the parameter to analyze
param=[0.2 log(2)/5 0.5 log(2)/80 1 1 5 5]; % reaction parameters 
param2=[0.2 log(2)/5 0.5 log(2)/20]; % reaction parameters
trsl=[0.15,0.25,0.40,0.55,0.70,0.90]; % translation rates to analyze

X=[0 1 1 0 0 0 0 0;
   0 0 0 1 0 0 1 1;
   1 0 0 0 1 0 1 0;
   0 0 0 0 0 1 0 1]; % Variables of each reaction.
X2=[1 0 0 0;
    0 1 1 0;
    0 0 0 1]; % Variables of each reaction.
S=[1 -1 0 0 0 0 0 0;
   0 0 1 -1 1 1 -1 -1;
   0 0 0 0 1 0 -1 0;
   0 0 0 0 0 1 0 -1]; % Stoichiometry of each reaction.
S2=[0 0 0 0;
    1 -1 0 0;
    0 0 1 -1]; % Stoichiometry of each reaction.
neg=[0 0 0 0 1 1 0 0]; % 1 if variable's coefficient is negative in transition rate wi(x).
neg2=[0 0 0 0]; % 1 if variable's coefficient is negative in transition rate wi(x). 
Const=[0 0 0 0 50 50 0 0]; % constants for transition rates wi(x)'s.
Const2=[0 0 0 0]; % constants for transition rates wi(x)'s.
Norm=[10 20 10 10]; % Scaling factor for each variable (repressor)
Norm2=[10 10 20]; % Scaling factor for each variable (reporter)

%% Select molecules to analyze
AnaMol1=2; % repressor part
AnaMol=4; % always 4
AnaMol2=3; % reporter part

%% Preparation
AveTetR=zeros(6,2);
VarTetR=zeros(6,2);
CVTetR=zeros(6,2);
AveGFP=zeros(6,2);
VarGFP=zeros(6,2);
CVGFP=zeros(6,2);

for i=1:6
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

    %% Solve SDP for TetR
        if SolveAve==1
            [avemax,avemin]=OrdervsAveSDP(A,H,mu,n,AnaMol1,Arow,Norm);
            AveTetR(i,1)=avemax
            AveTetR(i,2)=avemin
        end
        if SolveVar==1
            Varmax=OrdervsVarmaxSDP(A,H,mu,n,AnaMol1,Arow,Norm);
            VarTetR(i,1)=Varmax
            Varmin=OrdervsVarminSDP(A,H,mu,n,AnaMol1,Arow,Norm,avemax);
            VarTetR(i,2)=Varmin
        end
        if SolveCV==1
            CVmax=OrdervsCVmaxSDP(A,H,mu,n,AnaMol1,Arow,Norm);
            CVTetR(i,1)=CVmax
            CVmin=OrdervsCVminSDP(A,H,mu,n,AnaMol1,Arow,Norm,avemax);
            CVTetR(i,2)=CVmin
        end

    %% Solve SDP for DNA of GFP side
    if Analysis==1
        if SolveAve==1
            [avemax,avemin]=OrdervsAveSDP(A,H,mu,n,AnaMol,Arow,Norm);
        end
        if SolveVar+SolveCV>=1
            [M2ndmax,M2ndmin]=Ordervs2ndMomentSDP(A,H,mu,n,AnaMol,Arow,Norm);
        end

        %% Location of moments

        Acolumn2=momentlocation(2,n2);
        Arow2=[zeros(n2,1),momentlocation(mu+1,n2)];

        %% Determine A matrix

        J=[j1 j2 j3 j4];
        A2=Amatrix(2,n2,J,param2,X2,S2,neg2,Const2,Acolumn2,Arow2);
        A2=ANorm(A2,Norm2,Acolumn2,Arow2,2,n2); % Normalization

        %% Solve SDP
        if Analysis==1
            if SolveAve==1
                AveGFP(i,1)=avemax*(param2(1,1)*param2(1,3))/(param2(1,2)*param2(1,4))
                AveGFP(i,2)=avemin*(param2(1,1)*param2(1,3))/(param2(1,2)*param2(1,4))
            end
            if SolveVar==1
                varmax=OrdervsVarmaxSDP(A2,2,n2,AnaMol2,Arow2,Norm2,avemax,avemin,M2ndmax,M2ndmin);
                VarGFP(i,1)=varmax
                varmin=OrdervsVarminSDPGFP(A2,2,n2,AnaMol2,Arow2,Norm2,avemax,avemin,M2ndmax,M2ndmin,param2);
                VarGFP(i,2)=varmin
            end
            if SolveCV==1
                CVmax=OrdervsCVmaxSDPGFP(A2,2,n2,AnaMol2,Arow2,Norm2,avemax,avemin,M2ndmax,M2ndmin);
                CVGFP(i,1)=CVmax
                CVmin=OrdervsCVminSDPGFP(A2,2,n2,AnaMol2,Arow2,Norm2,avemax,avemin,M2ndmax,M2ndmin,param2);
                CVGFP(i,2)=CVmin
            end
        end              



    end
end
clear H
save SteadyStateAnalysis_RepressorReporter_Relation.mat