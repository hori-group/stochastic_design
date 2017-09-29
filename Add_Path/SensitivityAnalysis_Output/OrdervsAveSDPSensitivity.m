function [avemax,avemin]=OrdervsAveSDPSensitivity(A,As,H,mu,n,AnaMol,Arow,Norm,Upparam,Lwparam)

    At=zeros((3*size(H,1))^2,2*nchoosek(mu+n+1,mu+1)-1);
    H0=[H zeros(size(H,1),2*size(H,1),nchoosek(mu+n+1,mu+1))];
    Hu=[zeros(size(H,1),size(H,1),nchoosek(mu+n+1,mu+1)) Upparam*H zeros(size(H,1),size(H,1),nchoosek(mu+n+1,mu+1))];
    Hl=[zeros(size(H,1),2*size(H,1),nchoosek(mu+n+1,mu+1)) -Lwparam*H];
    Hs=[zeros(size(H,1),size(H,1),nchoosek(mu+n+1,mu+1)) zeros(size(H,1),size(H,1),nchoosek(mu+n+1,mu+1)) zeros(size(H,1),size(H,1),nchoosek(mu+n+1,mu+1));
        zeros(size(H,1),size(H,1),nchoosek(mu+n+1,mu+1)) -H zeros(size(H,1),size(H,1),nchoosek(mu+n+1,mu+1));
        zeros(size(H,1),size(H,1),nchoosek(mu+n+1,mu+1)) zeros(size(H,1),size(H,1),nchoosek(mu+n+1,mu+1)) H];
    H=[H0;Hu;Hl];
    
    for i=1:nchoosek(mu+n+1,mu+1)-1
        At(:,i)=-vec(H(:,:,i+1));
    end
    for i=1:nchoosek(mu+n+1,mu+1)
        At(:,nchoosek(mu+n+1,mu+1)+i-1)=-vec(Hs(:,:,i));
    end
    c=zeros(2*nchoosek(mu+n+1,mu+1)-1,1);
    Xk=zeros(n,1);
    Xk(AnaMol,1)=1;
    for i=1:nchoosek(mu+n+1,mu+1)-1
        if Arow(:,i+1)==Xk
            c(i,1)=1;
            location=i;
        end
    end
    u=A(:,1);
    A(:,1)=[];
    A=[A As];
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
    avemin=y(location,1)*Norm(1,AnaMol); % Minimum

    btt=c;
    [x,y,info]=sedumi(Att,btt,ctt,K,pars);
    info
    avemax=y(location,1)*Norm(1,AnaMol); % Maximum 

end