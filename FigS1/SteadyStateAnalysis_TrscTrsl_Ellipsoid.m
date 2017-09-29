hold on
%% Solve Mode   On:1 Off:0
Analysis=1;

%% define parameters
n=2; % The number of molecular species
mu=10; % The truncation order of the system.
i1=1; % The number of wi=ki reactions
i2=3; % The number of wi=kix reactions
i3=0; % The number of wi=kix(x-1) reactions
i4=0; % The number of wi=kix1x2 reactions
param=[0.2 log(2)/5 0.2 log(2)/20];%

X=[0 1 1 0;
   0 0 0 1]; % Choose variables
S=[1 -1 0 0;
   0 0 1 -1];
neg=[0 0 0 0];
Const=[20 0 0 0];
Norm=[30 200];

%% Select molecules to analyze
AnaMol=[1,2];

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
    [avemaxmRNA,aveminmRNA]=OrdervsAveSDP(A,H,mu,n,AnaMol(1,1),Arow,Norm);
    [avemaxProtein,aveminProtein]=OrdervsAveSDP(A,H,mu,n,AnaMol(1,2),Arow,Norm);
    varmaxmRNA=OrdervsVarmaxSDP(A,H,mu,n,AnaMol(1,1),Arow,Norm);
    varmaxProtein=OrdervsVarmaxSDP(A,H,mu,n,AnaMol(1,2),Arow,Norm);
    [DMax,CovMatrix]=DetMax(A,H,mu,n,AnaMol,Arow,Norm)
    [V,D]=eig(CovMatrix);
    LargeEig=D(2,2)
    SmallEig=D(1,1)
    LargeEigVec=V(:,2)
    SmallEigVec=V(:,1)
    MainAxis=sqrt(LargeEig);
    SubAxis=sqrt(SmallEig);
end

%% Draw ellipsoid
sin=LargeEigVec(2)/sqrt(LargeEigVec(1)^2+LargeEigVec(2)^2);
cos=LargeEigVec(1)/sqrt(LargeEigVec(1)^2+LargeEigVec(2)^2);
t=0:1:360;
a=MainAxis;
b=SubAxis;
x1=[a*cosd(t);
b*sind(t)];
R=[cos -sin;
sin cos];
X=R*x1;
figure(1)
plot(X(1,:)+(avemaxmRNA+aveminmRNA)/2,X(2,:)+(avemaxProtein+aveminProtein)/2,'-')