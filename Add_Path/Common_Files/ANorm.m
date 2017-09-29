function A=ANorm(A,Norm,Acolumn,Arow,mu,n)
    for i=1:nchoosek(mu+n,mu)-1 % Column count
        for j=1:nchoosek(mu+n+1,mu+1) % Row count
            for k=1:n
                A(i,j)=A(i,j)*Norm(1,k)^(Arow(k,j)-Acolumn(k,i));
            end
        end
    end
end