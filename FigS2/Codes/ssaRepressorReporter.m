clear;
cel=10000;
tEnd=1440;
x0 = [0;0;50;50;0;0];
k=[0.2 log(2)/5 0.4 log(2)/80 5 1 0.2 log(2)/5 0.5 log(2)/20];

S=[1 -1 0 0 0 0 0 0 0 0 0 0;
   0 0 1 -1 -1 1 -1 1 0 0 0 0;
   0 0 0 0 -1 1 0 0 0 0 0 0;
   0 0 0 0 0 0 -1 1 0 0 0 0;
   0 0 0 0 0 0 0 0 1 -1 0 0;
   0 0 0 0 0 0 0 0 0 0 1 -1]; %stoichiometry matrix (size of S = (# of molecules) x (# of reactions))
w = inline('[p(1)*x(3);p(2)*x(1);p(3)*x(1);p(4)*x(2);p(5)*x(2)*x(3);p(6)*(50-x(3));p(5)*x(2)*x(4);p(6)*(50-x(4));p(7)*x(4);p(8)*x(5);p(9)*x(5);p(10)*x(6)]','x','p');  % propensity vector (size of w = (# of reactions))

rand('state',sum(100*clock));

%% Run Gillespie
X=zeros(6,cel);
parfor i=1:cel
    [t,x] = RepressorReporternrm(S,w,tEnd,x0,k);
    kmax=size(x,2);
    X(:,i)=x(:,kmax);
end
RepressorMean=mean(X(2,:))
RepressorCV=sqrt(var(X(2,:)))/mean(X(2,:))
ReporterCV=sqrt(var(X(6,:)))/mean(X(6,:))

h=histogram(X(2,:))
h.Normalization = 'probability'
h.BinWidth = 2
xlabel('Copy number of repressor')
ylabel('Probability')

save ssaRepressorReporter.mat
