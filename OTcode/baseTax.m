addpath C:/Users/sebastian/Mosek/10.2/toolbox/r2017a
load segerpcAPclass ; dim = 50; %dimension of data
numPat = length(stride); %number of samples
celltypes= unique(cellnames);
numcells = length(celltypes);

%we create a list of pairwise sample combinations where the first sample
%is the reference sample. The first sample will always have the max number
%of clusters
allpats= 1:numPat;
refPats=find(stride==max(stride)); %Indices for samples with the most cell types.
refcombs= []; 
for i=1:length(refPats)
    pats = allpats(~ismember(allpats, refPats(i)));
    refcombs=[refcombs; combvec(refPats(i),pats)'];
end
%Note that the first column in refcombs is for the reference and the 
%second one for the unassigned sample. Notice also that the reference
%sample always moves in ascending order. the unassigned sample also moves
%in ascending order while the skipping the current reference sample.
lambda2=0;
lambda = .075;
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------
cellClusters = [];
perfMat = [];
CCMat=[];
predLabels=[];
cperf=[];
%Y=Y-1; needed if Y is made up of 1s and 2s
for i=1:length(refcombs)
    nclust1=stride(refcombs(i,1));  %number of clusters for individual 1
    nclust2=stride(refcombs(i,2));
    start1 = sum(stride(1:refcombs(i,1)-1))+1;   %how many columns to ignore +1 
    start2 = sum(stride(1:refcombs(i,2)-1))+1;
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
    gammaijcol =  gammaij./ max(abs(gammaij), [], 1);
    gammaijrow =  gammaij./ max(abs(gammaij), [], 2);
    gammaij = (gammaijcol+gammaijrow)/2;
    ground1=string(cellnames(start1:start1+nclust1-1));
    ground2=string(cellnames(start2:start2+nclust2-1));
    
    [~, Indices] = max(gammaij);  %row indices
    [~, assignments] = ismember(ground2, ground1);
    cellClusters = [cellClusters, Indices];
    if mod(i,numPat-1)==0 %applies when all samples have been aligned with current reference sample 
        cell_names = string(cellnames); %this copy will be modified so that clusters from ref are delted
        cells_pat=cellspat;
        newstride=stride;
        newY=Y;
        index1 =sum(stride(1:refcombs(i,1)-1));
        cell_names(index1+1:index1+max(stride))= [];
        cells_pat(refcombs(i,1))= [];
        newstride(refcombs(i,1)) = [];
        newY(refcombs(i,1)) = [];
        %cell_names= replace(cell_names,setdiff(unique(string(cellnames)), ground1),"unassigned");
        failcount = 0;
        fail = 0;
        for j=1:max(stride)
            PLindices = find(cellClusters==j); %indices for predicted labels
            predictedLabels = string(cell_names(PLindices)); % for the clusters 
            %of cluster-of-clusters i, we have different pred labels.
            %predLabel=mode(categorical(predictedLabels));
            predLabel = ground1(j); %predicted cell type for cluster of clusters j
            predLabels=[predLabels,predLabel];
            failcount = failcount+ sum(predictedLabels~=string(predLabel));
            propLabels = ww(PLindices);
            for k=1:length(PLindices)
                if predictedLabels(k)~=string(predLabel)
                    cellsp=cells_pat(min(find(cumsum(newstride)>= PLindices(k)))); %this is the 
                    %total number of cells in the patient where cell type j is from. 
                    fail=fail+round(propLabels(k)*cellsp);
                end
            end
    
        end
        [RI, ARI] = randindex(string(cellClusters),string(cell_names));
        perfMat = [perfMat;[ARI,1-failcount/length(cell_names),1-double(fail)/double(sum(cells_pat))]];        
        CCMat=[CCMat;cellClusters];
        %{d 
        %CLASSIFICATION
        predCellMat = zeros(numPat-1,nclust1);
        for j = 1:(numPat-1)
            start = sum(newstride(1:j-1))+1;
            predCellMat(j,1:newstride(j))=cellClusters(start:start+newstride(j)-1); 
        end
        propsmatNew = zeros(numPat-1, nclust1);
        for m=1:(numPat-1)
            start = sum(newstride(1:m-1))+1;
            props=ww(start:start+newstride(m)-1);
            patclusters=predCellMat(m,:); 
            propindex = 0; %
            for n=1:nclust1
                if ismember(patclusters(n), 1:nclust1) ==1
                    propindex = propindex+1;
                    if propsmatNew(m,patclusters(n))==0
                        propsmatNew(m,patclusters(n)) = props(propindex);
                    else
                        propsmatNew(m,patclusters(n)) = propsmatNew(m,patclusters(n))+props(propindex);
                    end
                end
            end
        end
        [accu, auc]=L1outRF(newY,numPat-1, propsmatNew,200);
        cperf = [cperf;[accu,auc]];
        %}
        cellClusters = [];
    end
end
perfMatOTRMC=perfMat;
OTRMCmean=mean(perfMat);
OTRMCrf = mean(cperf)
%OTRMCmean=mean(rates);

%####################################################################################
%####################################################################################3
rates=zeros(length(refcombs),1);
cellClusters = [];
perfMat = [];
CCMat=[];
predLabels=[];
cperf=[];
%Y=Y-1;
for i=1:length(refcombs)
    nclust1=stride(refcombs(i,1));  %number of clusters for individual 1
    nclust2=stride(refcombs(i,2));
    start1 = sum(stride(1:refcombs(i,1)-1))+1;   %how many columns to ignore +1 
    start2 = sum(stride(1:refcombs(i,2)-1))+1;
    ms1 = supp(1:dim,start1:start1+nclust1-1);
    ms2 = supp(1:dim,start2:start2+nclust2-1);
    vars1=supp(dim+1:dim+dim^2,start1:start1+nclust1-1);
    vars2 = supp(dim+1:dim+dim^2,start2:start2+nclust2-1);
    p1 = ww(start1:start1+nclust1-1);
    p2 = ww(start2:start2+nclust2-1);
    cost = real(CostMat(ms1,ms2,vars1,vars2,nclust1,nclust2));
    cost = cost/max(cost,[],"all");
    [~,res]=OT(cost,p1,p2);
    xx=res.sol.itr.xx;
    gammaij=reshape(xx,[nclust1,nclust2]);
    gammaijcol =  gammaij./ max(abs(gammaij), [], 1);
    gammaijrow =  gammaij./ max(abs(gammaij), [], 2);
    gammaij = (gammaijcol+gammaijrow)/2;
    ground1=string(cellnames(start1:start1+nclust1-1));
    ground2=string(cellnames(start2:start2+nclust2-1));
    
    [~, Indices] = max(gammaij);  %row indices
    [~, assignments] = ismember(ground2, ground1);
    cellClusters = [cellClusters, Indices];
    if mod(i,numPat-1)==0
        cell_names = string(cellnames); %this copy will be modified so that clusters from ref are delted
        cells_pat=cellspat;
        newstride=stride;
        newY=Y;
        index1 =sum(stride(1:refcombs(i,1)-1));
        cell_names(index1+1:index1+max(stride))= [];
        cells_pat(refcombs(i,1))= [];
        newstride(refcombs(i,1)) = [];
        newY(refcombs(i,1)) = [];
        %cell_names= replace(cell_names,setdiff(unique(string(cellnames)), ground1),"unassigned");
        failcount = 0;
        fail = 0;
        for j=1:max(stride)
            PLindices = find(cellClusters==j); %indices for predicted labels
            predictedLabels = string(cell_names(PLindices)); % for the clusters 
            %of cluster-of-clusters i, we have different pred labels.
            %predLabel=mode(categorical(predictedLabels));
            predLabel = ground1(j); %predicted cell type for cluster of clusters j
            predLabels=[predLabels,predLabel];
            failcount = failcount+ sum(predictedLabels~=string(predLabel));
            propLabels = ww(PLindices);
            for k=1:length(PLindices)
                if predictedLabels(k)~=string(predLabel)
                    cellsp=cells_pat(min(find(cumsum(newstride)>= PLindices(k)))); %this is the 
                    %total number of cells in the patient where cell type j is from. 
                    fail=fail+round(propLabels(k)*cellsp);
                end
            end
    
        end
        [RI, ARI] = randindex(string(cellClusters),string(cell_names));
        perfMat = [perfMat;[ARI,1-failcount/length(cell_names),1-double(fail)/double(sum(cells_pat))]];        
        CCMat=[CCMat;cellClusters];
        %{ d
        %CLASSIFICATION
        predCellMat = zeros(numPat-1,nclust1);
        for j = 1:(numPat-1)
            start = sum(newstride(1:j-1))+1;
            predCellMat(j,1:newstride(j))=cellClusters(start:start+newstride(j)-1); 
        end
        propsmatNew = zeros(numPat-1, nclust1);
        for m=1:(numPat-1)
            start = sum(newstride(1:m-1))+1;
            props=ww(start:start+newstride(m)-1);
            patclusters=predCellMat(m,:); 
            propindex = 0; %
            for n=1:nclust1
                if ismember(patclusters(n), 1:nclust1) ==1
                    propindex = propindex+1;
                    if propsmatNew(m,patclusters(n))==0
                        propsmatNew(m,patclusters(n)) = props(propindex);
                    else
                        propsmatNew(m,patclusters(n)) = propsmatNew(m,patclusters(n))+props(propindex);
                    end
                end
            end
        end
        [accu, auc]=L1outRF(newY,numPat-1, propsmatNew,200);
        cperf = [cperf;[accu,auc]];
        %}
        cellClusters = [];
    end
end
OTmean=mean(perfMat)
OTrf = mean(cperf)








