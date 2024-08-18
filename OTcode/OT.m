function [r,res]=OT(cost,p1,p2)
clear prob
p1=p1/sum(p1);
p2=p2/sum(p2);
m1=size(cost,1); %number of clusters in the first partition
m2=size(cost,2); %number of clusters in the second partition

coststack=reshape(cost,m1*m2,1); 
prob.c=coststack;
m1m2=length(coststack);

%Variable bounds
prob.blx=zeros(m1m2,1); 
prob.bux=[];

%Constraints.
%row sum to p1
a=zeros(m1+m2,m1m2);
idx=1:m2:m1m2;
for i=1:m1
    a(i,idx(i):idx(i)+m2-1)=a(i,idx(i):idx(i)+m2-1)+1;
end
%column sum to p2
for j=m1+1:m1+m2
    idx=j-m1:m1:m1m2;
    a(j,idx)=a(j,idx)+1;
end
prob.a=sparse(a);

prob.blc=[transpose(p1);transpose(p2)];
prob.buc=[transpose(p1);transpose(p2)];
%param.MSK_IPAR_INFEAS_REPORT_AUTO = 'MSK_ON';
param.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; 
%the interior point optimizer is used
param.MSK_IPAR_LOG = 1;
[r,res]=mosekopt('minimize',prob,param);
end