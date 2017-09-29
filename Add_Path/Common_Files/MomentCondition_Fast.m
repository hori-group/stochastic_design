function H=MomentCondition_Fast(mu,n,Arow)

%% The case of even order

    if mod(mu,2)==0 
        
    %% Preparation
        H=zeros((n+1)*nchoosek(mu/2+n,mu/2),(n+1)*nchoosek(mu/2+n,mu/2),nchoosek(mu+n+1,mu+1));
        H(1,1,1)=1;
        Xbase=Arow;
        Xbase(:,nchoosek(mu/2+n,mu/2)+1:nchoosek(mu+n+1,mu+1))=[];
        X=zeros(n,n+1);
        for k=1:n
            if k>0
                X(k,k+1)=1;
            end
        end
        
    %% Locate indices
        for local=0:n  
%            local
            for i=1:nchoosek(mu/2+n,mu/2)
                for j=1:nchoosek(mu/2+n,mu/2)
                    Hindex=Xbase(:,i)+Xbase(:,j)+X(:,local+1);
                    InSum=sum(Hindex);
                    if InSum>=1
                        for l=nchoosek(InSum+n-1,InSum-1)+1:nchoosek(InSum+n,InSum)
                            if Hindex==Arow(:,l)
                                H(local*nchoosek(mu/2+n,mu/2)+i,local*nchoosek(mu/2+n,mu/2)+j,l)=1;
                            end
                        end
                    end
                end
            end
        end        
     end

%% The case of odd order

    if mod(mu,2)==1 
        
    %% Preparation
        H=zeros(n*nchoosek((mu-1)/2+n,(mu-1)/2)+nchoosek((mu+1)/2+n,(mu+1)/2),n*nchoosek((mu-1)/2+n,(mu-1)/2)+nchoosek((mu+1)/2+n,(mu+1)/2),nchoosek(mu+n+1,mu+1));
        H(1,1,1)=1;
        Xbase=Arow;
        Xbase(:,nchoosek((mu+1)/2+n,(mu+1)/2)+1:nchoosek(mu+n+1,mu+1))=[];
        X=zeros(n,n+1);
        for k=1:n
            if k>0
                X(k,k+1)=1;
            end
        end
        
    %% Locate indices
        for local=0:n
            if local==0
                for i=1:nchoosek((mu+1)/2+n,(mu+1)/2)
                    for j=1:nchoosek((mu+1)/2+n,(mu+1)/2)
                        Hindex=Xbase(:,i)+Xbase(:,j)+X(:,local+1);
                        InSum=sum(Hindex);
                        if InSum>=1
                            for l=nchoosek(InSum+n-1,InSum-1)+1:nchoosek(InSum+n,InSum)
                                if Hindex==Arow(:,l)
                                    H(i,j,l)=1;
                                end
                            end
                        end
                    end
                end
            else
                for i=1:nchoosek((mu-1)/2+n,(mu-1)/2)
                    for j=1:nchoosek((mu-1)/2+n,(mu-1)/2)
                        Hindex=Xbase(:,i)+Xbase(:,j)+X(:,local+1);
                        InSum=sum(Hindex);
                        if InSum>=1
                            for l=nchoosek(InSum+n-1,InSum-1)+1:nchoosek(InSum+n,InSum)
                                if Hindex==Arow(:,l)
                                    H((local-1)*nchoosek((mu-1)/2+n,(mu-1)/2)+nchoosek((mu+1)/2+n,(mu+1)/2)+i,(local-1)*nchoosek((mu-1)/2+n,(mu-1)/2)+nchoosek((mu+1)/2+n,(mu+1)/2)+j,l)=1;
                                end
                            end
                         end
                    end
                end
            end
        end
    end
end