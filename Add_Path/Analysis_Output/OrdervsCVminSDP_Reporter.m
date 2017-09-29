function CVmin=OrdervsCVminSDP_Reporter(A,H,mu,n,AnaMol,Arow,Norm,avemax)

    At=zeros(size(H,1)^2,nchoosek(mu+n+1,mu+1)-1);
    for i=1:nchoosek(mu+n+1,mu+1)-1
        At(:,i)=-vec(H(:,:,i+1));
    end
    c=zeros(nchoosek(mu+n+1,mu+1)-1,1);
    Xk=zeros(n,1);
    Xk(AnaMol,1)=1;
    for i=1:nchoosek(mu+n+1,mu+1)-1
        if Arow(:,i+1)==2*Xk
            c(i,1)=1;
            location=i;
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
    clear A;
    clear H;
    pars.eps=10^(-9);
    pars.theta=0.01
    pars.beta=0.9
    pars.stepdif=1
    
    pars.cg.refine=10;
    pars.cg.maxiter=50;
    pars.cg.stagtol=5*10^(-20);
    pars.cg.restol=5*10^(-10);
    
    pars.chol.canceltol=10^(-20);
    
    pars.chol.maxuden=50;
    [x,y,info]=sedumi(Att,btt,ctt,K,pars);
    info
    CVmin2=y(location,1)*Norm(1,AnaMol)^2/avemax^2-1; % Minimum^2
    if CVmin2<0
        CVmin=0;
    else
        CVmin=sqrt(CVmin2);
    end
    
end