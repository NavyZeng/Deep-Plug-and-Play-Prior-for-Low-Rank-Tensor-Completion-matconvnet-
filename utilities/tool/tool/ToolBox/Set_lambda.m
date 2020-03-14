function lambda = Set_lambda(Nway,alpha) 
N = length(Nway);
lambda = 0;
for i=1:N-1
        for j=i+1:N
            d = prod(Nway)/(Nway(i)*Nway(j));
            temp=alpha(i,j)/sqrt(max(Nway(i),Nway(j))*d);
            lambda = lambda + temp;
        end
end

