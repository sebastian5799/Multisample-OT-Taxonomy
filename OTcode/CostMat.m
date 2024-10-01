
function cost = CostMat(ms1,ms2,vars1,vars2,nclust1,nclust2)
    cost = zeros(nclust1,nclust2);
    d=length(ms1(:,1)); %dimension of data
    for i = 1:nclust1
        m1 = ms1(:,i);  %mean vector for some cluster of pt 1.
        s1 = reshape(vars1(:,i),[d,d]);
        for j = 1:nclust2
            m2 = ms2(:,j);
            s2 = reshape(vars2(:,j),[d,d]);
            v= m1-m2;
            d2= v.'*v+ trace(s1+s2-(2*sqrtm(sqrtm(s1)*s2*sqrtm(s1)))); %Wasserstein distance
            cost(i,j) = cost(i,j)+sqrt(d2);
        end
    end 
end

