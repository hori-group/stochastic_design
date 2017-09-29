clear
%% Solve Mode   On:1 Off:0
Analysis=1;
    SolveAve=1; % Compute the upper and lower bounds of Average.
    SolveVar=0; % Compute the upper and lower bounds of Variance. Note that SolveAve must be also 1 if SolveVar is 1.
    SolveCV=1; % Compute the upper and lower bounds of CV. Note that SolveAve must be also 1 if SolveCV is 1.

%% Define parameters
n=3; % The number of variables (molecular species)
i1=0; % The number of wi=ki reactions (zero-th order)
i2=5; % The number of wi=kix reactions (first order)
i3=0; % The number of wi=kix(x-1) reactions (second order)
i4=1; % The number of wi=kix1x2 reactions (second order)
MaxOrder=10; % The maximum truncation order to finish analysis (note that order starts with 2nd)

param=[0.2 log(2)/5 0.5 log(2)/20 1 5]; % reaction parameters

X=[0 1 1 0 0 0;
   0 0 0 1 0 1;
   1 0 0 0 1 1]; % Variables of each reaction. 
S=[1 -1 0 0 0 0;
   0 0 1 -1 1 -1;
   0 0 0 0 1 -1]; % Stoichiometry of each reaction.
neg=[0 0 0 0 1 0]; % 1 if variable's coefficient is negative in transition rate wi(x). 
Const=[0 0 0 0 50 0]; % Constant term of transition rate wi(x). 
Norm=[5 20 5]; % Scaling factor for each variable

%% Choose Molecule to analyze
AnaMol=2;

%% Preparation (Memory Allocation)
Ave=zeros(MaxOrder-1,2);
Var=zeros(MaxOrder-1,2);
CV=zeros(MaxOrder-1,2);
Order=zeros(MaxOrder-1,1);

%% Calculation

for i=1:MaxOrder-1
    mu=i+1
    Order(i,1)=mu;
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
        if SolveVar==1
            varmax=OrdervsVarmaxSDP(A,H,mu,n,AnaMol,Arow,Norm);
            Var(i,1)=varmax
            varmin=OrdervsVarminSDP(A,H,mu,n,AnaMol,Arow,Norm,avemax);
            Var(i,2)=varmin
        end
        if SolveCV==1
            CVmax=OrdervsCVmaxSDP(A,H,mu,n,AnaMol,Arow,Norm);
            CV(i,1)=CVmax
            CVmin=OrdervsCVminSDP(A,H,mu,n,AnaMol,Arow,Norm,avemax);
            CV(i,2)=CVmin
        end
    end
end

%% Outputs

figure(1)
plot(Order,Ave(:,2),'--o',Order,Ave(:,1),'--x');
legend('Minimum Value','Maximum Value')
xlabel('Truncation order','FontName','Calibri','FontSize',20)
ylabel('Mean copy number of repressor','FontName','Calibri','FontSize',20)
set(gca,'FontName','Calibri','FontSize',20); 
figure(2)
plot(Order,Var(:,2),'--o',Order,Var(:,1),'--x');
legend('Minimum Value','Maximum Value')
xlabel('Truncation order','FontName','Calibri','FontSize',20)
ylabel('Variance of repressor','FontName','Calibri','FontSize',20)
set(gca,'FontName','Calibri','FontSize',20); 
figure(3)
plot(Order,CV(:,2),'--o',Order,CV(:,1),'--x');
legend('Minimum Value','Maximum Value')
xlabel('Truncation order','FontName','Calibri','FontSize',20)
ylabel('CV of repressor','FontName','Calibri','FontSize',20)
set(gca,'FontName','Calibri','FontSize',20); 


save SteadyStateAnalysis_OrdervsStatistics.mat