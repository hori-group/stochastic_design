function CVmax=OrdervsCVmaxSDP(A,H,mu,n,AnaMol,Arow,Norm)

    A=[A zeros(nchoosek(mu+n,mu)-1,1)];
    u=A(:,1);
    A(:,1)=[];
    At=zeros((size(H,1)+2)^2,nchoosek(mu+n+1,mu+1));
    c=zeros(nchoosek(mu+n+1,mu+1),1);
    c(nchoosek(mu+n+1,mu+1),1)=1;
    Xk=zeros(n,1);
    Xk(AnaMol,1)=1;
    for i=1:nchoosek(mu+n+1,mu+1)
        if Arow(:,i)==Xk
            location=i;
        end
    end
    for i=1:nchoosek(mu+n+1,mu+1)
        if Arow(:,i)==2*Xk
            location2=i;
        end
    end
    H=horzcat(H,zeros(size(H,1),2,nchoosek(mu+n+1,mu+1)));
    H=vertcat(H,zeros(2,size(H,2),nchoosek(mu+n+1,mu+1)));
    H(size(H,1)-1,size(H,1)-1,location2)=1;
    H(size(H,1)-1,size(H,1),location)=1;
    H(size(H,1),size(H,1)-1,location)=1;
    HV=zeros(size(H,1),size(H,1));
    HV(size(H,1),size(H,1))=1;
    for i=1:nchoosek(mu+n+1,mu+1)-1
        At(:,i)=-vec(H(:,:,i+1));
    end
    At(:,nchoosek(mu+n+1,mu+1))=-vec(HV);
    Att=[-A;At];
    btt=-c;
    ct=vec(H(:,:,1));
    ctt=[u;ct];
    K.f=size(A,1);
    K.s=size(H,1);
    pars.eps=0;
    pars.theta=0.01
    pars.beta=0.9
    pars.stepdif=1
    
    pars.cg.maxiter=500;
    pars.cg.refine=10;
    pars.cg.stagtol=5*10^(-20);
    pars.cg.restol=5*10^(-10);
    
    pars.chol.canceltol=10^(-20);
    
    pars.chol.maxuden=2000;
    
    [x,y,info]=sedumi(Att,btt,ctt,K,pars);
    info
    gammamin=y(nchoosek(mu+n+1,mu+1),1); % gamma(min)
    CVmax=(1/gammamin-1)^(1/2);

end