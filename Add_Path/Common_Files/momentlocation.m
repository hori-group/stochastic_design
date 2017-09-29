function muCount=momentlocation(mu,n)
    Count=zeros(n,(mu+1)^n-1);
    for i=1:(mu+1)^n-1
        Rownum=i;
        for j=1:n
            M=mod(Rownum,(mu+1)^(n-j));
            Count(n-j+1,i)=(Rownum-M)/((mu+1)^(n-j));
            if Count(n-j+1,i)>=1
                Rownum=Rownum-Count(n-j+1,i)*(mu+1)^(n-j);
            end
        end
    end
    muCount=zeros(n,nchoosek(mu+n,mu)-1);
    c=1;
    for i=1:mu
        for j=1:(mu+1)^n-1
            if sum(Count(:,j))==i
                muCount(:,c)=Count(:,j);
                c=c+1;
            end
        end
    end
end
