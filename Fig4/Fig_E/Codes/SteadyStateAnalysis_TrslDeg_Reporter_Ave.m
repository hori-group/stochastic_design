clear
%% Solve Mode   On:1 Off:0
Analysis=1;
    SolveAve=1; % Compute the upper and lower bounds of Average.

%% define parameters
n=4; % The number of molecular species in the "repressor" part of the circuit.
n2=3; % The number of molecular species in the "reporter" part of the circuit.
mu=7; % Truncation order of the moments.
i1=0; % The number of wi=ki reactions in the "repressor" part of the circuit.
i2=6; % The number of wi=kix reactions in the "repressor" part of the circuit.
i3=0; % The number of wi=kix(x-1) reactions in the "repressor" part of the circuit.
i4=2; % The number of wi=kix1x2 reactions in the "repressor" part of the circuit.
j1=0; % The number of wi=ki reactions in the "reporter" part of the circuit.
j2=4; % The number of wi=kix reactions in the "reporter" part of the circuit.
j3=0; % The number of wi=kix(x-1) reactions in the "reporter" part of the circuit.
j4=0; % The number of wi=kix1x2 reactions in the "reporter" part of the circuit.
Num=20; % The number of plots for each parameter
param=[0.2 log(2)/5 0.5 log(2)/20 1 1 5 5]; % reaction parameters for the "repressor" part.
param2=[0.2 log(2)/5 0.5 log(2)/20]; % reaction parameters for the "reporter" part.
trsl=0.05:0.05:0.05+0.05*(Num-1); % translation rates to analyze
halfT=5:5:5+5*(Num-1); % Half-life time of protein to analyze
de=zeros(1,Num); 
for i=1:Num
    de(1,i)=log(2)/halfT(1,i); % Transform half-life times to degradation rates
end

X=[0 1 1 0 0 0 0 0;
   0 0 0 1 0 0 1 1;
   1 0 0 0 1 0 1 0;
   0 0 0 0 0 1 0 1]; % Variables of each reaction (repressor part)
X2=[1 0 0 0;
    0 1 1 0;
    0 0 0 1]; % Variables of each reaction (reporter part)
S=[1 -1 0 0 0 0 0 0;
   0 0 1 -1 1 1 -1 -1;
   0 0 0 0 1 0 -1 0;
   0 0 0 0 0 1 0 -1]; % Stoichiometry of each reaction (repressor part)
S2=[0 0 0 0;
    1 -1 0 0;
    0 0 1 -1]; % Stoichiometry of each reaction (reporter part)
neg=[0 0 0 0 1 1 0 0]; % 1 if variable's coefficient is negative in transition rate wi(x). (repressor part)
neg2=[0 0 0 0]; % 1 if variable's coefficient is negative in transition rate wi(x). (repressor part)
Const=[0 0 0 0 50 50 0 0]; % Constants for transition rate wi(x)'s. (repressor part)
Const2=[0 0 0 0]; % Constants for transition rate wi(x)'s (reporter part)
Norm=[10 20 10 10]; % Scaling factor for each variable (repressor part)
Norm2=[10 10 20]; % Scaling factor for each variable (reporter part)

%% Select molecules to analyze
AnaMol=4; % reporessor part
AnaMol2=3; % reporter part

%% Memory allocation
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
                    Ave(i,j,1)=avemax*(param2(1,1)*param2(1,3))/(param2(1,2)*param2(1,4))
                    Ave(i,j,2)=avemin*(param2(1,1)*param2(1,3))/(param2(1,2)*param2(1,4))
                end
            end    
        end
    end
end
 
figure(1)
pcolor(trsl,halfT,Ave(:,:,1));
xlabel('Translation Gain')
ylabel('Half Degradation Time')
zlabel('Average')
    
save SteadyStateAnalysis_TrslDeg_Reporter_Ave.mat