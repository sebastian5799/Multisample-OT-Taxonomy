function[ARI, clustA, cellA] = TaxFunction(dataFile,dim,lambda)
load(dataFile)
numPat = length(stride); %number of samples
celltypes= unique(cellnames);
numcells = length(celltypes);
%We will now apply the OTrmc algorithm to each patient with herself.
lambda2=0;
matchSelf = cell(numPat,1);
for i =1:numPat
    nclust=stride(i);
    start = sum(stride(1:i-1))+1;
    ms = supp(1:dim,start:start+nclust-1);
    vars = supp(dim+1:dim^2+dim,start:start+nclust-1);
    p = ww(start:start+nclust-1);
    cost = CostMat(ms,ms,vars,vars,nclust,nclust);
    cost = cost/max(cost,[],"all");
    [~,res]=otrmcl1(cost,lambda,lambda2,p,p);
    xx=res.sol.itr.xx;
    gammaij=reshape(xx(nclust+nclust+1:nclust+nclust+nclust*nclust),[nclust,nclust]);
    %normalization
    gammaijcol =  gammaij./ max(abs(gammaij), [], 1);
    gammaijrow =  gammaij./ max(abs(gammaij), [], 2);
    matchSelf{i} = (gammaijcol+gammaijrow)/2;
end

%We will compute the matching matrices for the different sample
%combinations
combs = nchoosek(1:numPat,2);
matchM = cell(length(combs),1);
for i=1:length(combs)
    nclust1=stride(combs(i,1));  %number of clusters for individual 1
    nclust2=stride(combs(i,2));
    start1 = sum(stride(1:combs(i,1)-1))+1;   %how many columns to ignore +1 
    start2 = sum(stride(1:combs(i,2)-1))+1;
    ms1 = supp(1:dim,start1:start1+nclust1-1);
    ms2 = supp(1:dim,start2:start2+nclust2-1);
    vars1=supp(dim+1:dim+dim^2,start1:start1+nclust1-1);
    vars2 = supp(dim+1:dim+dim^2,start2:start2+nclust2-1);
    p1 = ww(start1:start1+nclust1-1);
    p2 = ww(start2:start2+nclust2-1);
    cost = CostMat(ms1,ms2,vars1,vars2,nclust1,nclust2);
    cost = cost/max(cost,[],"all");
    [~,res]=otrmcl1(cost,lambda,lambda2,p1,p2);
    xx=res.sol.itr.xx;
    gammaij=reshape(xx(nclust1+nclust2+1:nclust1+nclust2+nclust1*nclust2),[nclust1,nclust2]);
    %normalization
    gammaijcol =  gammaij./ max(abs(gammaij), [], 1);
    gammaijrow =  gammaij./ max(abs(gammaij), [], 2);
    matchM{i} = (gammaijcol+gammaijrow)/2;
end

%WE will now combine all of the different matching matrices into one big
%matching matrix. 
Biggammaij = blkdiag(matchSelf{:}); % this functions takes a finite  list of matrices 
                                    % and creates a larger matrix witht he
                                    % matrices as diagonal blocks. {:}
                                    % spits out all the matrices from the
                                    % cell array
startPat=1;
endPat = numPat-1;
for j =1:(numPat-1)
    nclust1=stride(j);  %number of clusters for individual j
    start1 = sum(stride(1:j-1))+1;   %how many columns to ignore +1
    endClust=start1+nclust1-1;
    Biggammaij(start1:endClust,endClust+1:end)=Biggammaij(start1:endClust,endClust+1:end)+ [matchM{startPat:endPat}];
    Biggammaij(endClust+1:end,start1:endClust)=Biggammaij(endClust+1:end,start1:endClust)+ [matchM{startPat:endPat}]';
    startPat = endPat+1;
    endPat = startPat+numPat-2-j;
end

%Now we will create a distance matrix
A = real(-log(Biggammaij));
A= A/ max(max(A));
B = squareform(A);
tree=linkage(B,'ward'); 
cellClusters=cluster(tree,"MaxClust", numcells);
%dendrogram(tree,'ColorThreshold',24)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now we will compute the Accuracy based on cell (as opposed to cluster). To
%do this we need to know the number of cells that are present in each
%cluster. Another option is to use the proportion of each cluster which is
%stored in the ww variable. Along the way we can also compute an cluster 
%accuracy again with a different approach. 
failcount = 0;
fail = 0;
for i=1:numcells
    PLindices = find(cellClusters==i); %indices for predicted labels
    predictedLabels = string(cellnames(PLindices)); % for the clusters 
    %of cluster-of-clusters i, we have different pred labels.
    predLabel=mode(categorical(predictedLabels));
    failcount = failcount+ sum(predictedLabels~=string(predLabel));
    propLabels = ww(PLindices);
    cumstride = cumsum(stride);
    for j=1:length(PLindices)
        if predictedLabels(j)~=string(predLabel)
            cellsp=cellspat(min(find(cumsum(stride)>= PLindices(j)))); %this is the 
            %total number of cells in the patient where cell type j is from. 
            fail=fail+round(propLabels(j)*cellsp);
        end
    end
    
end
%CLuster accuracy
clustA=1-failcount/length(cellnames);
%Cell accuracy
cellA=1-double(fail)/double(sum(cellspat));
[RI, ARI] = randindex(string(cellClusters),string(cellnames));
end