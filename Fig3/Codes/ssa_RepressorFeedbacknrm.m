function [tReaction,xTraj] = ssa_RepressorFeedbacknrm(S,w,tEnd,x0,k)
%
% S: stoichio metry (#mol. species \times #reaction)
% prop: propensities (#reactions)
% tEnd: endtime
% x0: initial value (#mol. species)
%
p = k;
t = 0;
x = x0;
tReaction = [0];
xTraj = x;

while t <= tEnd
   wx = w(x,p);
   sumw = sum(wx);
   t = t + 1/sumw*log(1/rand);
   if t > tEnd
       break;
   end
   tReaction = [tReaction, t];
   
   ref = sumw*rand;
   comp = wx(1);

   for i=1:length(wx)
       if comp >= ref
           x = x + S(:,i);
           xTraj = [xTraj, x];
           break;
       end
       comp = comp + wx(i+1);
   end
end
end
