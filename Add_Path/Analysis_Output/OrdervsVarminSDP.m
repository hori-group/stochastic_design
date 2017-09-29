function varmin=OrdervsVarminSDP(A,H,mu,n,AnaMol,Arow,Norm,avemax)

    At=zeros(size(H,1)^2,nchoosek(mu+n+1,mu+1)-1);
    for i=1:nchoosek(mu+n+1,mu+1)-1
        At(:,i)=-vec(H(:,:,i+1));
    end
    c=zeros(nchoosek(mu+n+1,mu+1)-1,1);
    Xk=zeros(n,1);
    Xk(AnaMol,1)=1;
    for i=1:nchoosek(mu+n+1,mu+1)-1
        if Arow(:,i+1)==Xk
            c(i,1)=-avemax/Norm(1,AnaMol);
            location1=i;
        end
    end
    for i=1:nchoosek(mu+n+1,mu+1)-1
        if Arow(:,i+1)==2*Xk
            c(i,1)=1;
            location2=i;
        end
    end
    u=A(:,1);
    A(:,1)=[];
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
    
    pars.chol.maxuden=4000;
    [x,y,info]=sedumi(Att,btt,ctt,K,pars);
    info
    varmin=y(location2,1)*Norm(1,AnaMol)^2-avemax*y(location1,1)*Norm(1,AnaMol); % Minimum
    if varmin<0
        varmin=0;
    end
    
end