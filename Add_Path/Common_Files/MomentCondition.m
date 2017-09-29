function H=MomentCondition(mu,n,Arow)

%% The case of even order

    if mod(mu,2)==0 
        
    %% Preparation
        H=zeros((n+1)*nchoosek(mu/2+n,mu/2),(n+1)*nchoosek(mu/2+n,mu/2),nchoosek(mu+n+1,mu+1));
        Hindex=zeros((n+1)*nchoosek(mu/2+n,mu/2),(n+1)*nchoosek(mu/2+n,mu/2),n);
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
            for k=1:n
                for i=1:nchoosek(mu/2+n,mu/2)
                    for j=1:nchoosek(mu/2+n,mu/2)
                        Hindex(local*nchoosek(mu/2+n,mu/2)+i,local*nchoosek(mu/2+n,mu/2)+j,k)=Xbase(k,i)+Xbase(k,j)+X(k,local+1);
                    end
                end
            end
        end
        
    %% Preparation
        Arow3=zeros(1,nchoosek(mu+n+1,mu+1),n);
        for i=1:nchoosek(mu+n+1,mu+1)
            for j=1:n
                Arow3(1,i,j)=Arow(j,i);
            end
        end
        
    %% Put moments
        H(1,1,1)=1;
        for k=2:nchoosek(mu+n+1,mu+1)
            for i=1:(n+1)*nchoosek(mu/2+n,mu/2)
                for j=1:(n+1)*nchoosek(mu/2+n,mu/2)
                    if Hindex(i,j,:)==Arow3(1,k,:)
                        H(i,j,k)=1;
                    end
                end
            end
        end
        
    end

%% The case of odd order

    if mod(mu,2)==1 
        
    %% Preparation
        H=zeros(n*nchoosek((mu-1)/2+n,(mu-1)/2)+nchoosek((mu+1)/2+n,(mu+1)/2),n*nchoosek((mu-1)/2+n,(mu-1)/2)+nchoosek((mu+1)/2+n,(mu+1)/2),nchoosek(mu+n+1,mu+1));
        Hindex=zeros(n*nchoosek((mu-1)/2+n,(mu-1)/2)+nchoosek((mu+1)/2+n,(mu+1)/2),n*nchoosek((mu-1)/2+n,(mu-1)/2)+nchoosek((mu+1)/2+n,(mu+1)/2),n);
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
                for k=1:n
                    for i=1:nchoosek((mu+1)/2+n,(mu+1)/2)
                        for j=1:nchoosek((mu+1)/2+n,(mu+1)/2)
                            Hindex(i,j,k)=Xbase(k,i)+Xbase(k,j)+X(k,1);
                        end
                    end
                end
            else
                for k=1:n
                    for i=1:nchoosek((mu-1)/2+n,(mu-1)/2)
                        for j=1:nchoosek((mu-1)/2+n,(mu-1)/2)
                            Hindex((local-1)*nchoosek((mu-1)/2+n,(mu-1)/2)+nchoosek((mu+1)/2+n,(mu+1)/2)+i,(local-1)*nchoosek((mu-1)/2+n,(mu-1)/2)+nchoosek((mu+1)/2+n,(mu+1)/2)+j,k)=Xbase(k,i)+Xbase(k,j)+X(k,local+1);
                        end
                    end
                end
            end
        end
        
    %% Preparation
        Arow3=zeros(1,nchoosek(mu+n+1,mu+1),n);
        for i=1:nchoosek(mu+n+1,mu+1)
            for j=1:n
                Arow3(1,i,j)=Arow(j,i);
            end
        end
        
    %% Put moments
        H(1,1,1)=1;
        for k=2:nchoosek(mu+n+1,mu+1)
            for i=1:n*nchoosek((mu-1)/2+n,(mu-1)/2)+nchoosek((mu+1)/2+n,(mu+1)/2)
                for j=1:n*nchoosek((mu-1)/2+n,(mu-1)/2)+nchoosek((mu+1)/2+n,(mu+1)/2)
                    if Hindex(i,j,:)==Arow3(1,k,:)
                        H(i,j,k)=1;
                    end
                end
            end
        end
        
    end
end