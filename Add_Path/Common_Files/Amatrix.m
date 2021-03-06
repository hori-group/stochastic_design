function Aall=Amatrix(mu,n,I,param,X,S,neg,Const,Acolumn,Arow)

A=ones(nchoosek(mu+n,mu)-1,nchoosek(mu+n+1,mu+1),sum(I)); % Memory allocation

%% Reaction i1: wi=ki reactions
    if I(1,1)>0
        for reac=1:I(1,1)
            for i=1:nchoosek(mu+n,mu)-1 % Column count
                for j=1:nchoosek(mu+n+1,mu+1) % Row count
                    for k=1:n
                        if Acolumn(k,i)>0 && Acolumn(k,i)>=Arow(k,j)
                            A(i,j,reac)=A(i,j,reac)*nchoosek(Acolumn(k,i),Arow(k,j))*S(k,reac)^(Acolumn(k,i)-Arow(k,j));
                        elseif Acolumn(k,i)<Arow(k,j)
                            A(i,j,reac)=0;
                        end
                    end
                    if Acolumn(:,i)==Arow(:,j)
                        A(i,j,reac)=0;
                    end
                    A(i,j,reac)=param(1,reac)*Const(1,reac)*A(i,j,reac);
                end
            end
        end
    end
    Isum=I(1,1);
    
%% Reaction i2: wi=kix reactions
    if I(1,2)>0
        for reac=1+Isum:I(1,2)+Isum
            A2=ones(nchoosek(mu+n,mu)-1,nchoosek(mu+n+1,mu+1)); % Memory allocation
            for i=1:nchoosek(mu+n,mu)-1 % Column count
                for j=1:nchoosek(mu+n+1,mu+1) % Row count
                    for k=1:n
                        if Acolumn(k,i)>0 && Acolumn(k,i)+X(k,reac)>=Arow(k,j) && Arow(k,j)>=X(k,reac)
                            A(i,j,reac)=A(i,j,reac)*nchoosek(Acolumn(k,i),Arow(k,j)-X(k,reac))*S(k,reac)^(Acolumn(k,i)-Arow(k,j)+X(k,reac));
                        elseif Acolumn(k,i)+X(k,reac)<Arow(k,j) || Arow(k,j)<X(k,reac)
                            A(i,j,reac)=0;
                        end
                    end
                    if Acolumn(:,i)+X(:,reac)==Arow(:,j)
                        A(i,j,reac)=0;
                    end
                    A(i,j,reac)=param(1,reac)*A(i,j,reac);
                    if neg(1,reac)==1
                        for k=1:n
                            if Acolumn(k,i)>0 && Acolumn(k,i)>=Arow(k,j)
                                A2(i,j)=A2(i,j)*nchoosek(Acolumn(k,i),Arow(k,j))*S(k,reac)^(Acolumn(k,i)-Arow(k,j));
                            elseif Acolumn(k,i)<Arow(k,j)
                                A2(i,j)=0;
                            end
                        end
                        if Acolumn(:,i)==Arow(:,j)
                            A2(i,j)=0;
                        end
                        A2(i,j)=param(1,reac)*Const(1,reac)*A2(i,j);
                    end
                end
            end
            if neg(1,reac)==1
                A(:,:,reac)=A2-A(:,:,reac);
            end
        end
    end
    Isum=I(1,2)+Isum;
    
%% Reaction i3: wi=kix(x-1) reactions
    if I(1,3)>0
        for reac=1+Isum:I(1,3)+Isum
            A2=ones(nchoosek(mu+n,mu)-1,nchoosek(mu+n+1,mu+1)); % Memory allocation
            for i=1:nchoosek(mu+n,mu)-1 % Column count
                for j=1:nchoosek(mu+n+1,mu+1) % Row count
                    for k=1:n
                        if Acolumn(k,i)>0 && Acolumn(k,i)+X(k,reac)*2>=Arow(k,j) && Arow(k,j)>=X(k,reac)*2
                            A(i,j,reac)=A(i,j,reac)*nchoosek(Acolumn(k,i),Arow(k,j)-X(k,reac)*2)*S(k,reac)^(Acolumn(k,i)-Arow(k,j)+X(k,reac)*2);
                        elseif Acolumn(k,i)+X(k,reac)*2<Arow(k,j) || Arow(k,j)<X(k,reac)*2
                            A(i,j,reac)=0;
                        end
                    end
                    if Acolumn(:,i)+X(:,reac)*2==Arow(:,j)
                        A(i,j,reac)=0;
                    end
                    A(i,j,reac)=param(1,reac)*A(i,j,reac);
                    for k=1:n
                        if Acolumn(k,i)>0 && Acolumn(k,i)+X(k,reac)>=Arow(k,j) && Arow(k,j)>=X(k,reac)
                            A2(i,j)=A2(i,j)*nchoosek(Acolumn(k,i),Arow(k,j)-X(k,reac))*S(k,reac)^(Acolumn(k,i)-Arow(k,j)+X(k,reac));
                        elseif Acolumn(k,i)+X(k,reac)<Arow(k,j) || Arow(k,j)<X(k,reac)
                            A2(i,j)=0;
                        end
                    end
                    if Acolumn(:,i)+X(:,reac)==Arow(:,j)
                        A2(i,j)=0;
                    end
                    A2(i,j)=param(1,reac)*A2(i,j);
                end
            end
            A(:,:,reac)=A(:,:,reac)-A2;
        end
    end
    Isum=I(1,3)+Isum;


%% Reaction i4: wi=kix1x2 reactions
    if I(1,4)>0
        for reac=1+Isum:I(1,4)+Isum
            for i=1:nchoosek(mu+n,mu)-1 % Column count
                for j=1:nchoosek(mu+n+1,mu+1) % Row count
                    for k=1:n
                        if Acolumn(k,i)>0 && Acolumn(k,i)+X(k,reac)>=Arow(k,j) && Arow(k,j)>=X(k,reac)
                            A(i,j,reac)=A(i,j,reac)*nchoosek(Acolumn(k,i),Arow(k,j)-X(k,reac))*S(k,reac)^(Acolumn(k,i)-Arow(k,j)+X(k,reac));
                        elseif Acolumn(k,i)+X(k,reac)<Arow(k,j) || Arow(k,j)<X(k,reac)
                            A(i,j,reac)=0;
                        end
                    end
                    if Acolumn(:,i)+X(:,reac)==Arow(:,j)
                        A(i,j,reac)=0;
                    end
                    A(i,j,reac)=param(1,reac)*A(i,j,reac);
                end
            end
        end
    end
%% Decide A matrix
    Aall=sum(A,3);
end

