function [r,res]=otrmcl1(cost,lambda,lambda2,p1,p2)

%cost[m1,m2] is the cost distance between any pair of clusters
%p1[m1] is the marginal for the first distribution, column
%vector p2[m2] is the marginal for the second distribution,
%column vector.

  
%  cost=[0.1,0.2,1.0;0.8,0.8,0.1]; p1=[0.4;0.6];
%  p2=[0.2;0.3;0.5];lambda=1.0;

clear prob;

m1=size(cost,1); %number of clusters in the first partition
m2=size(cost,2); %number of clusters in the second partition

%Variables: 1 to m1: gap variables for the marginal 
%constraints of distribution 1; m1+1 to m1+m2: gap variables
%for the marginal constraints of distribution 2; m1+m2+1 to
%m1+m2+m1*m2: matching weights gamma_ij

coststack=reshape(cost,m1*m2,1); %column-wise scan
prob.c=[lambda*ones(m1+m2,1);coststack];
%prob.c=[1000*ones(m1,1);lambda*ones(m2,1);coststack];  %why?

%same dimension
prob.qosubi=1:(m1+m2);
prob.qosubj=1:(m1+m2);
prob.qoval=ones(m1+m2,1)*2*lambda2; 

%same dimension
prob.blx=zeros(m1+m2+m1*m2,1); 
prob.bux=inf*ones(m1+m2+m1*m2,1);

%same number of rows
prob.a=zeros(2*m1+2*m2+1,m1+m2+m1*m2);
prob.blc=zeros(2*m1+2*m2+1,1);
prob.buc=zeros(2*m1+2*m2+1,1);


omegasub=1:(m1+m2);
gammasub=reshape(m1+m2+1:m1+m2+m1*m2,[m1,m2]);

for i=1:m1
    prob.a(i,omegasub(i))=-1.0;
    prob.a(i+m1,omegasub(i))=1.0;
    for j=1:m2
        prob.a(i,gammasub(i,j))=1.0;
        prob.a(i+m1, gammasub(i,j))=1.0;

        prob.buc(i)=p1(i);
        prob.blc(i)=-inf;

        prob.buc(i+m1)=inf;
        prob.blc(i+m1)=p1(i);

    end
end


for j=1:m2
	prob.a(j+2*m1,omegasub(m1+j))=-1.0;
    prob.a(j+2*m1+m2,omegasub(m1+j))=1.0;
    for i=1:m1
        prob.a(j+2*m1,gammasub(i,j))=1.0;
        prob.a(j+2*m1+m2, gammasub(i,j))=1.0;
    end

    prob.buc(2*m1+j)=p2(j);
    prob.blc(2*m1+j)=-inf;

    prob.buc(2*m1+j+m2)=inf;
    prob.blc(2*m1+j+m2)=p2(j);

end

for i=1:m1
	for j=1:m2
		prob.a(2*m1+2*m2+1,gammasub(i,j))=1.0;
    end
end
prob.blc(2*m1+2*m2+1)=1.0;
prob.buc(2*m1+2*m2+1)=1.0;

param.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; 
%the interior point optimizer is used
param.MSK_IPAR_LOG = 0;
       
[r,res]=mosekopt('minimize',prob,param);



