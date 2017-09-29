clear
%% Solve Mode   On:1 Off:0
Analysis=1;
    SolveAve=1; % Compute the upper and lower bounds of Average.
    SolveVar=0; % Compute the upper and lower bounds of Variance. Note that SolveAve must be also 1 if SolveVar is 1.
    SolveCV=1; % Compute the upper and lower bounds of CV. Note that SolveAve must be also 1 if SolveCV is 1.

%% define parameters
n=4; % The number of molecular species in the repressor part of the circuit
n2=3; % The number of molecular species in the reporter part of the circuit
mu=6; % truncation order of the moments.
i1=1; % The number of wi=ki reactions (repressor part)
i2=6; % The number of wi=kix reactions 
i3=0; % The number of wi=kix(x-1) reactions
i4=2; % The number of wi=kix1x2 reactions 
j1=0; % The number of wi=ki reactions (reporter part)
j2=4; % The number of wi=kix reactions
j3=0; % The number of wi=kix(x-1) reactions
j4=0; % The number of wi=kix1x2 reactions

% Note that 1st and 6th reactions are actually the same reactions but separately defined to analyze the effect of perturbations to DNA total copy number.
param=[1 0.2 log(2)/5 0.4 log(2)/80 1 1 5 5]; % parameters
param2=[0.2 log(2)/5 0.5 log(2)/20]; % parameters
PerRate=0.05; % magnitude of perturbation (0.05 = 5% of the nominal value)
X=[0 0 1 1 0 0 0 0 0;
   0 0 0 0 1 0 0 1 1;
   0 1 0 0 0 1 0 1 0;
   0 0 0 0 0 0 1 0 1]; % Variables of each reaction (repressor)
X2=[1 0 0 0;
    0 1 1 0;
    0 0 0 1]; % Variables of each reaction (reporter)
S=[0 1 -1 0 0 0 0 0 0;
   1 0 0 1 -1 1 1 -1 -1;
   1 0 0 0 0 1 0 -1 0;
   0 0 0 0 0 0 1 0 -1]; % Stoichiometry of each reaction (repressor)
S2=[0 0 0 0;
    1 -1 0 0;
    0 0 1 -1]; % Stoichiometry of each reaction (reporter)
neg=[0 0 0 0 0 1 1 0 0]; % 1 if variable's coefficient is negative in transition rate wi(x).
neg2=[0 0 0 0]; % 1 if variable's coefficient is negative in transition rate wi(x). 
Const=[50 0 0 0 0 0 50 0 0]; % constants for transition rates wi(x)'s
Const2=[0 0 0 0]; % constants for transition rates wi(x)'s
Norm=[20 50 10 10]; % Scaling factor for each variable (repressor)
Norm2=[10 20 50]; % Scaling factor for each variable (reporter)

%% Select molecules to analyze
AnaMol=[2,4];
AnaMol2=3;

%% memory allocation
paramSAv=1:1:5;
Ave=zeros(2,length(paramSAv),2);
Var=zeros(2,length(paramSAv),2);
CV=zeros(2,length(paramSAv),2);
ave0max=zeros(2,1);
ave0min=zeros(2,1);
var0max=zeros(2,1);
CV0max=zeros(2,1);


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
       [ave0max(1,1),ave0min(1,1)]=OrdervsAveSDP(A,H,mu,n,AnaMol(1,1),Arow,Norm)
    end
    if SolveVar==1
        var0max(1,1)=OrdervsVarmaxSDP(A,H,mu,n,AnaMol(1,1),Arow,Norm)
    end
    if SolveCV==1
        CV0max(1,1)=OrdervsCVmaxSDP(A,H,mu,n,AnaMol(1,1),Arow,Norm)
    end
end

j=1;
    for i=1:length(paramSAv)

        %% The upper and lower bounds of the parameter to analyze
        paramSA=i;
        Upparam=param(1,paramSA)*(1+PerRate)
        Lwparam=param(1,paramSA)*(1-PerRate)
        if i==3 || i==5
            Upparam=param(1,paramSA)*(1/(1-PerRate))
            Lwparam=param(1,paramSA)*(1/(1+PerRate))
        end
        paramA=param;
        paramA(1,paramSA)=0;
        noparam=zeros(1,i1+i2+i3+i4);
        noparam(1,paramSA)=1;

        %% Location of moments

        Acolumn=momentlocation(mu,n);
        Arow=[zeros(n,1),momentlocation(mu+1,n)];

        %% Determine A matrix

        I=[i1 i2 i3 i4];
        A=Amatrix(mu,n,I,paramA,X,S,neg,Const,Acolumn,Arow);
        A=ANorm(A,Norm,Acolumn,Arow,mu,n); % Normalization

        %% Determine As matrix

        I=[i1 i2 i3 i4];
        As=Amatrix(mu,n,I,noparam,X,S,neg,Const,Acolumn,Arow);
        As=ANorm(As,Norm,Acolumn,Arow,mu,n); % Normalization

        %% Stieltjes Moment Condition

        H=MomentCondition(mu,n,Arow);
            
        if j==1

            %% Solve SDP
            if Analysis==1
                if SolveAve==1
                   [avemax,avemin]=OrdervsAveSDPSensitivity(A,As,H,mu,n,AnaMol(1,j),Arow,Norm,Upparam,Lwparam);
                    Ave(1,i,j)=avemax
                    Ave(2,i,j)=avemin
                end
                if SolveVar==1
                    varmax=OrdervsVarmaxSDPSensitivity(A,As,H,mu,n,AnaMol(1,j),Arow,Norm,Upparam,Lwparam);
                    Var(1,i,j)=varmax
                end
                if SolveCV==1
                    CVmax=OrdervsCVmaxSDPSensitivity(A,As,H,mu,n,AnaMol(1,j),Arow,Norm,Upparam,Lwparam);
                    CV(1,i,j)=CVmax
                end
            end
        else
        end
        
    end
save SensitivityAnalysis.mat