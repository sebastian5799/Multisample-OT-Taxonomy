function [PerfMetrics, RFperf ] = RefAlign(stride, cellnames, ww, supp,cellspat, dim, method,Y)
    %This function is used to label the clusters of different samples using
    %the labeled clusters of one reference sample. We pick as the reference
    %sample that which has the largest number of clusters in it. In several
    %instances we may have multiple reference samples when several samples
    %have the same maximum number of clusters. If say we have 3 reference
    %samples, then we take the AVERAGE performance of the 3 labeling
    %procedures. This function offers 2 methods for labeling the clusters 
    %of samples: one is using OTRMC and the second is regular OT. 
    %Here is a review of the variables needed:
    %1)stride: vector storing the number of clusters for each sample
    %2)cellnames: string vector storing the ground truth cell labels for
    %  all clusters from all samples.
    %3)ww: vector storing the proportions of all clusters relative to the
    %  sample they come from. Same length as cellnames.
    %4)supp: matrix storing the mean vectors and covariance matrices of all
    %  clusters from all samples. 
    %5)cellspat: vector storing the number of cells in each sample. It has
    %  the same length as stride
    %6)dim: the dimension of the data, ie the length of the cluster mean
    %  vectors.
    %7)method: OTRMC vs OT
    %8)Y: vector storing sample labels, e.g. diseased vs healthy. Same
    %  length as stride. Used for Random Forest(RF) classification of
    %  samples
    %The outputs are two objects:
    %1)PerfMetrics: a vector containing the ARI (cluster ground truth 
    %  labels for all samples EXCEPT reference sample vs cluster predicted
    %  labels for all samples using ground truth reference sample labels), 
    %  Cluster Accuracy and Cell Accuracy.
    %2)RFperf: a vector storing SAMPLE (not cluster) classification accuracy and AUC.
    %---------------------------------------------------------------------
    RFperf = "Not Applicable"; %default RF output if Y is not provided.
    %---------------------------------------------------------------------
    %we create a list of pairwise sample combinations where the first sample
    %is the reference sample. 
    numPat = length(stride); %number of samples
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
    if nargin>7&&max(Y)==2 %Y is provided and is in 1,2 form
        Y=Y-1; %makes Y be in 0,1 form 
    end
    %----------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------
    cellClusters = [];
    perfMat = [];
    cperf=[];
    for i=1:length(refcombs)
        nclust1=stride(refcombs(i,1)); %number of clusters for sample 1
        nclust2=stride(refcombs(i,2));
        start1 = sum(stride(1:refcombs(i,1)-1))+1; %how many columns to ignore +1 
        start2 = sum(stride(1:refcombs(i,2)-1))+1;
        ms1 = supp(1:dim,start1:start1+nclust1-1); %mean vectors for clusters in sample 1
        ms2 = supp(1:dim,start2:start2+nclust2-1);
        vars1=supp(dim+1:dim+dim^2,start1:start1+nclust1-1); %vectorized vars for clusters in sample 1
        vars2 = supp(dim+1:dim+dim^2,start2:start2+nclust2-1);
        p1 = ww(start1:start1+nclust1-1); %cluster proportions in sample 1
        p2 = ww(start2:start2+nclust2-1);
        cost = CostMat(ms1,ms2,vars1,vars2,nclust1,nclust2);
        cost = cost/max(cost,[],"all"); 
        %%%%%%%%%
        if method == "OTRMC" %OTRMC used for reference alignment
            [~,res]=otrmcl1(cost,lambda,lambda2,p1,p2);
            xx=res.sol.itr.xx; %matching weight matrix
            gammaij=reshape(xx(nclust1+nclust2+1:nclust1+nclust2+nclust1*nclust2),[nclust1,nclust2]);
        else
            [~,res]=OT(cost,p1,p2);
            xx=res.sol.itr.xx;
            gammaij=reshape(xx,[nclust1,nclust2]);
        end
        %%%%%%%%%%
        gammaijcol =  gammaij./ max(abs(gammaij), [], 1);
        gammaijrow =  gammaij./ max(abs(gammaij), [], 2);
        gammaij = (gammaijcol+gammaijrow)/2;
        ground1=string(cellnames(start1:start1+nclust1-1));%Reference cell types  
        
        [~, Indices] = max(gammaij);  %row indices
        %[~, assignments] = ismember(ground2, ground1);
        cellClusters = [cellClusters, Indices]; %stores predictions for current reference sample
        if mod(i,numPat-1)==0 % i.e. All samples have been aligned with current reference sample 
            %we create copies of objects where the reference sample is
            %eliminated.
            cell_names = string(cellnames); 
            cells_pat=cellspat;
            newstride=stride;
            index1 =sum(stride(1:refcombs(i,1)-1)); %indices of reference clusters
            cell_names(index1+1:index1+max(stride))= [];
            cells_pat(refcombs(i,1))= [];
            newstride(refcombs(i,1)) = [];
            if nargin>7 %eliminate reference from sample labels
                newY=Y;
                newY(refcombs(i,1)) = [];
            end
            failcount = 0; %counts the number wrongly labeled clusters
            fail = 0; %counts the number of wrongly labeled cells
            for j=1:max(stride) %ie for each cluster in reference sample
                PLindices = find(cellClusters==j); %indices for clusters of metacluster j 
                predictedLabels = string(cell_names(PLindices)); % The clusters 
                %of metacluster j have different ground truth labels.
                predLabel = ground1(j); %predicted cell type for metacluster j
                failcount = failcount+ sum(predictedLabels~=string(predLabel));
                propLabels = ww(PLindices);
                for k=1:length(PLindices) %for each cluster k in metacluster j 
                    if predictedLabels(k)~=string(predLabel) %ie ground truth label of cluster 
                        % differs from that of the metacluster prediction 
                        cellsp=cells_pat(min(find(cumsum(newstride)>= PLindices(k)))); %this is the 
                        %total number of cells in the sample where cluster k is from. 
                        fail=fail+round(propLabels(k)*cellsp); %number of wrongly labeled cells
                    end
                end
        
            end
            %Performance Metrics: ARI, cluster labeling accuracy and
            %cell labeling accuracy.
            [~, ARI] = randindex(string(cellClusters),string(cell_names));
            perfMat = [perfMat;[ARI,1-failcount/length(cell_names),1-double(fail)/double(sum(cells_pat))]];        
            %CLASSIFICATION
            if nargin>7
                predCellMat = zeros(numPat-1,nclust1);
                for jj = 1:(numPat-1)
                    start = sum(newstride(1:jj-1))+1;
                    predCellMat(jj,1:newstride(jj))=cellClusters(start:start+newstride(jj)-1); 
                end
                propsmatNew = zeros(numPat-1, nclust1);
                for m=1:(numPat-1)
                    start = sum(newstride(1:m-1))+1;
                    props=ww(start:start+newstride(m)-1);
                    patclusters=predCellMat(m,:); 
                    propindex = 0; 
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
            end
            cellClusters = []; %we will use the next reference sample
        end
    end
    if isvector(perfMat) % if single reference present:
        PerfMetrics = perfMat;
    else  %several references present:
        PerfMetrics=mean(perfMat);
    end
    if nargin>7 %If sample classification was realized:
        RFperf = mean(cperf);
    end
end








