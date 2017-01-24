%Copy number variation analysis
%Import data
GSE57872_i = readtable('GSE57872_GBM_data_matrix.xlsx');
GSE57872_i = GSE57872_i(~all(isnan(GSE57872_i{:,3:end}),2),:);
chrom_minu = GSE57872_i{:,'chr'};
sortIDs = GSE57872_i.Properties.VariableNames(4:end);
exprgenesbycells = GSE57872_i{:,ismember(GSE57872_i.Properties.VariableNames,sortIDs)};
normalcells = {'MGH30_C12','MGH31_A08','MGH31_B07','MGH31_E04','MGH26_C03','MGH31_E09','MGH31_H07','MGH31_A04','MGH31_B05','MGH31_F10'};
geneschr = GSE57872_i.Var1;

geneexpr_i = exprgenesbycells(:,:);
% currently data is log2(TPM+1) - meanexpr
% would like to change to log2(TPM/10+1) - meanexpr

cellIDs = sortIDs(:,:);
tumorlabel = cell(size(cellIDs));
for i = 1:length(cellIDs);
    t = strsplit(cellIDs{i},'_');
    tumorlabel(i) = t(1);
end  
tumor = cellIDs(strcmp(tumorlabel,'MGH31'));


geneexpr = exprgenesbycells(:,ismember(cellIDs,tumor) | ismember(cellIDs,normalcells));
cellIDs = cellIDs(ismember(cellIDs,tumor) | ismember(cellIDs,normalcells));


chrom= chrom_minu(:,:);

% threshold relative gene expression values (log2(TPM+1))-meanofgene > 3 or < -3 
geneexpr_filt = geneexpr;
%first re-center genes to 0
for k = 1:size(geneexpr_filt,1);
    geneexpr_filt(k,:) = geneexpr_filt(k,:)-mean(geneexpr_filt(k,:));
end
geneexpr_filt(geneexpr_filt > 3) = 3;
geneexpr_filt(geneexpr_filt < -3) = -3;

%Moving average of 100 genes
CNV = zeros(size(geneexpr,1),size(geneexpr,2));
for k = 1:size(cellIDs,2)
    for i = 51:size(geneexpr,1)-50;
        CNV(i,k) = sum(geneexpr_filt((i-50):(i+50),k))/101;
    end
end

CNV = CNV(51:(size(geneexpr,1)-50),:);
chrom =chrom(51:(size(geneexpr,1)-50),:);
geneschr = geneschr(51:(size(geneexpr,1)-50),:);
% Center CNV vector of each cell at 0
for k = 1:size(cellIDs,2)
    CNV(:,k) = CNV(:,k) - median(CNV(:,k));
end

%clustergram(CNV,'Colormap',redbluecmap,'cluster',2)
%CNV_new = zeros(size(CNV,1),size(CNV,2)/2);
%for k = 1:size(CNV_new,2)
%    CNV_new(:,k) = CNV(:,k*2-1)-CNV(:,k*2);    
%end

%HeatMap(CNV_new,'Colormap',redbluecmap,'RowLabels',[1:size(chrom,1)],'RowLabelsColor',struct('labels',chrom,'color',{}))

%CNV_random = CNV(:,randsample(1:size(sortIDs,2),size(sortIDs,2)));
%CNV_randomminus=zeros(size(CNV,1),size(CNV,2)/2);
%for k = 1:size(CNV_new,2)
%    CNV_randomminus(:,k) = CNV_random(:,k*2-1)-CNV_random(:,k*2);    
%end

%HeatMap(CNV_randomminus,'Colormap',redbluecmap)


% subtract "normal" cluster
CNV_normal = mean(CNV(:,ismember(cellIDs,normalcells)),2);
CNV_norm = zeros(size(CNV));
for k = 1:size(cellIDs,2) %for each cell
    for i = 1:size(CNV,1)
        if CNV(i,k) > CNV_normal(i) + 0.3
            CNV_norm(i,k) = CNV(i,k) - CNV_normal(i);
        elseif CNV(i,k) < CNV_normal(i) -0.3
            CNV_norm(i,k) = CNV(i,k) - CNV_normal(i);
        else
            CNV_norm(i,k) = 0;
        end
    end
    
end
clear CNV_normal



C = clustergram(flipud(CNV_norm),'columnlabels',cellIDs,'cluster',2,'OptimalLeafOrder','true','Colormap',redbluecmap,'Linkage','ward','symmetric',true,'RowLabels',chrom,'DisplayRange',2.5);

%Classifying cells by subtype
genes = GSE57872_i.Var1;




%Classifying cells by subtype
%by the method of Broad
%first take average relative expression of each set of predictor genes minus
%average relative expression of all genes
markers_i = readtable('markerssubtype.txt','delimiter','\t','HeaderLines',1);
subtypes = markers_i.Properties.VariableNames;
subtype_score = zeros(size(geneexpr,2),size(subtypes,2));
submarkers = markers_i{:,:};
for t=1:size(submarkers,2) %for each subtype
    for j =1:size(geneexpr,2)  %for each cell
       subtype_score(j,t)= mean(geneexpr(ismember(genes,submarkers(:,t)),j)) - mean(geneexpr(:,j));
    end
end

subtype_IDs = zeros(size(geneexpr,2),size(subtypes,2));
subtype_robustness = zeros(size(geneexpr,2),size(subtypes,2));
for t=1:size(submarkers,2)
    %select 100 random gene sets
    randsets=cell(sum(ismember(genes,submarkers(:,t))),100);
    for i =1:100
        randsets(:,i)= randsample(genes,sum(ismember(genes,submarkers(:,t))));
    end
    %use rnadomsets to define 5% cutoff for expected subtype scores
    for j = 1:size(geneexpr,2)
        scoretemp=zeros(100,1);
        for i =1:100
            scoretemp(i) = mean(geneexpr(ismember(genes,randsets(:,i)),j)) - mean(geneexpr(:,j));
        end
        scoretemp=sort(scoretemp,'descend');
        if subtype_score(j,t) > scoretemp(5)
            subtype_IDs(j,t) = 1;
        elseif subtype_score(j,t) < scoretemp(95)
            subtype_IDs(j,t) = -1;
        end
        %Evaluate robustness of each subtype assignment
       % markertemp = submarkers(ismember(submarkers(:,t),genes),t);
       % for k=1:1000
       %     randset=datasample(markertemp,size(markertemp,1));
       %     if mean(geneexpr(ismember(genes,randset),j)) - mean(geneexpr(:,j))> scoretemp(5)
       %         subtype_robustness(j,t)=subtype_robustness(j,t)+1;
       %     end
       % end
    end
end



clusterorder = C.ColumnLabels;
resortI = zeros(size(clusterorder));%vector of IDs that I will sort subtype scores with
for i = 1:length(clusterorder)
    resortI(i) = find(strcmp(cellIDs,clusterorder(i)));
end
singlecellsubtypescores = subtype_score(resortI,:); % sort subtype scores to match clustering
singlecellsubtypeIDs = subtype_IDs(resortI,:);

H1= HeatMap(singlecellsubtypescores','Colormap',redbluecmap,'rowlabels',subtypes,'columnlabels',C.ColumnLabels,'DisplayRange',2);
cmap = [ 0.5725    0.7725    0.8706; 1 1 1; 0.8392    0.3765    0.3020];

H2= HeatMap(singlecellsubtypeIDs','Colormap',cmap,'rowlabels',subtypes,'columnlabels',C.ColumnLabels);

GroupIDs = zeros(size(clusterorder));
GroupIDs(ismember(clusterorder,C1.ColumnNodeNames)) = 1;
GroupIDs(ismember(clusterorder,C2.ColumnNodeNames)) = 2;
GroupIDs(ismember(clusterorder,C3.ColumnNodeNames)) = 3;
GroupIDs(ismember(clusterorder,C4.ColumnNodeNames)) = 4;



pvals = ones(size(singlecellsubtypescores,2),1);
for i =1:length(pvals)
    pvals(i) = kruskalwallis(singlecellsubtypescores(GroupIDs>0,i),GroupIDs(GroupIDs>0));
end

% as a control
% pull subtype genes' cnv scores
% compare distribution across cells

for s = 1:width(markers_i) %for each subtype
    I = ismember(geneschr,markers_i{:,s}); %find marker genes in CNV vector
    CNVtemp = CNV_norm(I,:); %pull out CNV vector for cells
    %test whether there is a difference between the average CNV for 
    
    
end


%for plotting where the chromosome boundries are
chromosome_vector = zeros(size(chrom));
for i = 1:length(chrom)
    if isequal(chrom{i},'X')
        chromosome_vector(i) = 23;
    else
        chromosome_vector(i) = str2double(chrom{i});
    end
end

H3 = HeatMap(flipud(chromosome_vector), 'symmetric','false','Colormap','colorcube');










