clear;
cel=10000; % Total number of cells (samples)
tEnd=1440; % Time to stop the simulation
x0 = [0;0]; % Initial state of the variables
k=[0.2*20 log(2)/5 0.2 log(2)/20]; % parameters

S=[1 -1 0 0;
   0 0 1 -1]; %stoichiometry matrix (size of S = (# of molecules) x (# of reactions))
w = inline('[p(1);p(2)*x(1);p(3)*x(1);p(4)*x(2)]','x','p');  % propensity vector (size of w = (# of reactions))

rand('state',sum(100*clock));

%% Run Gillespie
X=zeros(2,cel);
parfor i=1:cel
    [t,x] = ssaTrscTrslnrm(S,w,tEnd,x0,k);
    kmax=size(x,2);
    X(:,i)=x(:,kmax);
end

%% Statistics

meanmRNA=mean(X(1,:)) % Mean of mRNA copy number
meanProtein=mean(X(2,:)) % Mean of protein copy number
varmRNA=var(X(1,:)) % Variance of mRNA copy number
varProtein=var(X(2,:)) % Variance of protein copy number
SDmRNA=sqrt(varmRNA) % Standard deviation of mRNA copy number
SDProtein=sqrt(varProtein) % Standard deviation of protein copy number
CVmRNA=SDmRNA/meanmRNA % CV of mRNA copy number
CVProtein=SDProtein/meanProtein % CV of proteincopy number
fanomRNA=varmRNA/meanmRNA % Fano factor of mRNA copy number
fanoProtein=varProtein/meanProtein % Fano factor of protein copy number
covariance=cov(X') % Covariance of mRNA & protein copy number
correlation=covariance(1,2)/sqrt(varmRNA*varProtein) % correlation of mRNA & protein copy number

save mRNAvsProteinStatistics.mat


%% Figure Settings

MaxX1=50; 
MinX1=10; 
MaxX2=250;
MinX2=100;
IntX1=2; % size of each bin (x axis)
IntX2=10; % size of each bin (y axis)

%% Calculation for Output Figure

DataNumX1=(MaxX1-MinX1)/IntX1+1;
DataNumX2=(MaxX2-MinX2)/IntX2+1;
Prob=zeros(DataNumX1,DataNumX2);
mole1=MinX1:IntX1:MaxX1;
mole2=MinX2:IntX2:MaxX2;

for i=1:DataNumX1
    for j=1:DataNumX2
        count=0;
        for l=1:size(X,2)
            if X(1,l)>=MinX1+IntX1*(i-1) && MinX1+IntX1*i>X(1,l) && X(2,l)>=MinX2+IntX2*(j-1) && MinX2+IntX2*j>X(2,l)
                count=count+1;
            end
        end
        Prob(i,j)=count/cel;
    end
end

%% Output Figure

pcolor(mole1,mole2,Prob');
xlabel('mRNA copy number')
ylabel('protein copy number')
c=colorbar;
c.Label.String = 'Probability'