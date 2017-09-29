function [DMax,CovMatrix]=DetMax(A,H,mu,n,AnaMol,Arow,Norm)

    A=[A zeros(nchoosek(mu+n,mu)-1,4)];
    u=A(:,1);
    A(:,1)=[];
    At=zeros((size(H,1)+7)^2,nchoosek(mu+n+1,mu+1)+3);
    c=zeros(nchoosek(mu+n+1,mu+1)+3,1);
    c(nchoosek(mu+n+1,mu+1)+3,1)=-1;
    Xk1=zeros(n,1);
    Xk1(AnaMol(1),1)=1;
    Xk2=zeros(n,1);
    Xk2(AnaMol(2),1)=1;
    for i=1:nchoosek(mu+n+1,mu+1)
        if Arow(:,i)==Xk1
            location10=i;
        end
    end
    for i=1:nchoosek(mu+n+1,mu+1)
        if Arow(:,i)==2*Xk1
            location20=i;
        end
    end
    for i=1:nchoosek(mu+n+1,mu+1)
        if Arow(:,i)==Xk2
            location01=i;
        end
    end
    for i=1:nchoosek(mu+n+1,mu+1)
        if Arow(:,i)==2*Xk2
            location02=i;
        end
    end
    for i=1:nchoosek(mu+n+1,mu+1)
        if Arow(:,i)==Xk1+Xk2
            location11=i;
        end
    end
    H=horzcat(H,zeros(size(H,1),7,nchoosek(mu+n+1,mu+1)));
    H=vertcat(H,zeros(7,size(H,2),nchoosek(mu+n+1,mu+1)));
    
    H(size(H,1)-2,size(H,1)-2,1)=1;
    H(size(H,1)-3,size(H,1)-3,1)=1;
    H(size(H,1)-4,size(H,1)-4,1)=1;
    
    H(size(H,1)-2,size(H,1)-6,location10)=1;%*Norm(AnaMol(1));
    H(size(H,1)-6,size(H,1)-2,location10)=1;%*Norm(AnaMol(1));
    
    H(size(H,1)-6,size(H,1)-6,location20)=1;%*Norm(AnaMol(1))^2;
    
    H(size(H,1)-2,size(H,1)-5,location01)=1;%*Norm(AnaMol(2));
    H(size(H,1)-5,size(H,1)-2,location01)=1;%*Norm(AnaMol(2));
    
    H(size(H,1)-5,size(H,1)-5,location02)=1;%*Norm(AnaMol(2))^2;
    
    H(size(H,1)-6,size(H,1)-5,location11)=1;%*Norm(AnaMol(1))*Norm(AnaMol(2));
    H(size(H,1)-5,size(H,1)-6,location11)=1;%*Norm(AnaMol(1))*Norm(AnaMol(2));
    
    Hz11=zeros(size(H,1),size(H,1));
    Hz11(size(H,1)-1,size(H,1)-1)=1;
    Hz11(size(H,1)-4,size(H,1)-6)=1;
    Hz11(size(H,1)-6,size(H,1)-4)=1;
    
    Hz12=zeros(size(H,1),size(H,1));
    Hz12(size(H,1)-4,size(H,1)-5)=1;
    Hz12(size(H,1)-5,size(H,1)-4)=1;
    
    Hz22=zeros(size(H,1),size(H,1));
    Hz22(size(H,1),size(H,1))=1;
    Hz22(size(H,1)-3,size(H,1)-5)=1;
    Hz22(size(H,1)-5,size(H,1)-3)=1;
    
    HV=zeros(size(H,1),size(H,1));
    HV(size(H,1)-1,size(H,1))=1;
    HV(size(H,1),size(H,1)-1)=1;
    
    for i=1:nchoosek(mu+n+1,mu+1)-1
        At(:,i)=-vec(H(:,:,i+1));
    end
    At(:,nchoosek(mu+n+1,mu+1))=-vec(Hz11);
    At(:,nchoosek(mu+n+1,mu+1)+1)=-vec(Hz12);
    At(:,nchoosek(mu+n+1,mu+1)+2)=-vec(Hz22);
    At(:,nchoosek(mu+n+1,mu+1)+3)=-vec(HV);
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
    DMax=y(nchoosek(mu+n+1,mu+1)+3,1)^4*Norm(AnaMol(1))^2*Norm(AnaMol(2))^2;
    CovMatrix=[y(location20-1,1)-y(location10-1,1)^2,y(location11-1,1)-y(location10-1,1)*y(location01-1,1);
               y(location11-1,1)-y(location10-1,1)*y(location01-1,1),y(location02-1,1)-y(location01-1,1)^2];
    CovMatrix(1,1)=CovMatrix(1,1)*Norm(AnaMol(1))^2;
    CovMatrix(1,2)=CovMatrix(1,2)*Norm(AnaMol(1))*Norm(AnaMol(2));
    CovMatrix(2,1)=CovMatrix(2,1)*Norm(AnaMol(1))*Norm(AnaMol(2));
    CovMatrix(2,2)=CovMatrix(2,2)*Norm(AnaMol(2))^2;
    Dist=DMax-CovMatrix(1,1)*CovMatrix(2,2)+CovMatrix(2,1)*CovMatrix(1,2)
end