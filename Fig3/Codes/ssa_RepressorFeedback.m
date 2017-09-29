clear;
cel=10000; % Total number of cells (samples)
tEnd=1440; % Time to stop the simulation
x0 = [0;0;50]; % Initial state of the variables
k=[0.2 log(2)/5 0.5 log(2)/20 1 50 5]; % parameters
S=[1 -1 0 0 0 0;
   0 0 1 -1 1 -1;
   0 0 0 0 1 -1]; %stoichiometry matrix (size of S = (# of molecules) x (# of reactions))
w = inline('[p(1)*x(3);p(2)*x(1);p(3)*x(1);p(4)*x(2);p(5)*(p(6)-x(3));p(7)*x(2)*x(3)]','x','p');  % propensity vector (size of w = (# of reactions))

rand('state',sum(100*clock));

%% Run Gillespie
X=zeros(cel,1);
parfor i=1:cel
    [t,x] = ssa_RepressorFeedbacknrm(S,w,tEnd,x0,k);
    kmax=size(x,2);
    X(i,1)=x(2,kmax);
end
save ssa_RepressorFeedback.mat
h=histogram(X)
h.Normalization = 'probability'