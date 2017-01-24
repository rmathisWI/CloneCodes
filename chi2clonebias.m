%% Do clones have heterogeneous biases?
%compare to the median clone

% Load in data
reads_i = readtable('clonereads.txt','Delimiter','\t','ReadRowNames',1,'HeaderLines',1); %load in reads of clones
%sum postsort replicates
labeli = cellfun(@(x) {x(1:2)},reads_i.Properties.VariableNames);
label = unique(labeli);
reads = array2table(zeros(height(reads_i),length(label)),'RowNames',reads_i.Properties.RowNames,'VariableNames',label);
for t = 1:length(label)
    reads{:,label{t}} = sum(reads_i{:,strcmp(labeli,label{t})},2); %sum postsort replicates
    reads{isnan(reads{:,label{t}}),label{t}} = 0; %replace NaNs (not found) with 0s
end



%% convert reads to fractions
%divide by the library size
rnorm = zeros(size(reads_i));
for t = 1:width(reads_i);
    switch reads_i.Properties.VariableNames{t}(1)
        case 'E'
            Srt = 0.5968;
        case 'M'
            Srt = 0.4032;
    end    
    rnorm(:,t) = Srt.*(reads_i{:,t}./nansum(reads_i{:,t}));
end
rnorm(isnan(rnorm)) = 1E-6;
rnorm = array2table(rnorm,'VariableNames',reads_i.Properties.VariableNames,'RowNames',reads_i.Properties.RowNames);


%% calculate fraction epithelial
timepoints = cellfun(@str2num,(cellfun(@(x) {x(2)}, rnorm.Properties.VariableNames)));
states = cellfun(@(x) {x(1)}, rnorm.Properties.VariableNames);
frac_ep = zeros(height(rnorm),max(timepoints));
for t = 1:max(timepoints)
    E = sum(rnorm{:,timepoints==t & strcmp(states,'E')},2);
    C = sum(rnorm{:,timepoints==t},2);
    frac_ep(:,t) = E./C;
end

%% calculate log2(E/M)
L2EM_t = zeros(height(rnorm),max(timepoints));
timepoints = cellfun(@str2num,(cellfun(@(x) {x(2)}, rnorm.Properties.VariableNames)));
states = cellfun(@(x) {x(1)}, rnorm.Properties.VariableNames);
for t = 1:max(timepoints)
    E = sum(rnorm{:,timepoints==t & strcmp(states,'E')},2);
    M = sum(rnorm{:,timepoints==t & strcmp(states,'M')},2);
    L2EM_t(:,t) = log2(E./M);
end


%% Test if clones are sig. different from the population average

timepoints = cellfun(@str2num,(cellfun(@(x) {x(2)}, reads.Properties.VariableNames)));
states = cellfun(@(x) {x(1)}, reads.Properties.VariableNames);
expected = zeros(size(reads));
for t = 1:max(timepoints);
    E_obs = reads{:,timepoints==t & strcmp(states,'E')};
    M_obs = reads{:,timepoints==t & strcmp(states,'M')};
    obs = sum(reads{:,timepoints==t},2);
    %If same ratio as the parental population, proportion of reads for each
    %clone equals the ratio of library reads;
    r = sum(E_obs)/sum(obs);
    expected(:,timepoints==t & strcmp(states,'E')) = obs.*r;
    expected(:,timepoints==t & strcmp(states,'M')) = obs.*(1-r);
end


% determine chi squared
Xsqr = zeros(height(reads),max(timepoints));  
for t = 1:size(Xsqr,2);
    tp = timepoints == timepoints(t); %find columns of time points
    
   % for c = 1:size(Xsqr,1);
     %   Xsqr(c,t) = sum( (reads{c,tp} - expected(c,tp)).^2 ./ expected(c,tp));
   % end
    
    Xsqr(:,t) = sum(((reads{:,tp} - expected(:,tp)).^2)./ expected(:,tp),2); %sum of ((obs-exp)^2)/exp  
end

pvals_bias = chi2cdf(Xsqr,1,'upper');
s = size(pvals_bias);
FDR_bias = reshape(mafdr(pvals_bias(:),'BHFDR',1),s);
bias = all(FDR_bias < 0.05,2);
[~,i] = min(abs(frac_ep-.6));

%% test if reads are more than expected from contamination
%this requires computing the reads from contamination.
%calculate the probability a read from a cell of state S ends up in the
%not-S sorted population 
reads_i = readtable('clonereads.txt','Delimiter','\t','ReadRowNames',1,'HeaderLines',1); %load in reads of clones
labeli = cellfun(@(x) {x(1:2)},reads_i.Properties.VariableNames);
label = unique(labeli);
reads = array2table(zeros(height(reads_i),length(label)),'RowNames',reads_i.Properties.RowNames,'VariableNames',label);
for t = 1:length(label)
    reads{:,label{t}} = sum(reads_i{:,strcmp(labeli,label{t})},2); %sum postsort replicates
    reads{isnan(reads{:,label{t}}),label{t}} = 0; %replace NaNs (not found) with 0s
end

timepoints = cellfun(@str2num,(cellfun(@(x) {x(2)}, reads.Properties.VariableNames)));
states = cellfun(@(x) {x(1)}, reads.Properties.VariableNames);
%contamination_fraction = frac contam. in  [E1, E2, E3, M1, M2, M3];
cf_labels = {'E1','E2','E3','M1','M2','M3'};
contamination_fraction = [0.0488, 0.0263, 0.0260, 0.0648, 0.0441, 0.0353];
contam_reads = zeros(size(reads));
pvals_notcontam = ones(size(reads));
n=0;
for s = 1 : 2; % for each state
    for t = 1 : max(timepoints) %for each time point
        n=n+1;
        switch s
            case 1 % epithelial
                %Srt = 0.5968;
                Cor = 'E';  %correctly-sorted gate 
                InCor = 'M'; %incorrectly-sorted gate
            case 2 % mesenchymal
                %Srt = 0.4032;
                Cor = 'M';  %correctly-sorted gate 
                InCor = 'E'; %incorrectly-sorted gate
        end 
        
        cfr_Cor = contamination_fraction(strcmp(cf_labels,horzcat(Cor,num2str(t)))); % fraction contamination of sorted population of state S
        cfr_InCor = contamination_fraction(strcmp(cf_labels,horzcat(InCor,num2str(t)))); % fraction contamination of other sorted population
        i_Cor = strcmp(states,Cor) & timepoints == t; % ID of state S
        i_InCor = strcmp(states,InCor) & timepoints == t; % ID of other state
        
        S = sum(reads{:,i_Cor}) * (1-cfr_Cor); % reads for correctly sorted cells of state S
        S_prime = sum(reads{:,i_InCor}) * cfr_InCor; % reads for incorrectly sorted cells of state S
        pC = S_prime / (S + S_prime);  % proportion of reads from state S that are mis-sorted
        
        %now use binomial test to determine significance of more reads than expect.
        rds = sum(reads{:,timepoints == t},2);
        obs = reads{:,i_InCor};
        pvals_notcontam(:,n) = binocdf(obs,rds,pC,'upper');
    end
end

%make q values: turn pvals to vector, calculate Q values, return to matrix.
s = size(pvals_notcontam);
FDR = mafdr(pvals_notcontam(:),'BHFDR',1);
FDRvals_notcontam = reshape(FDR,s);

% how many clones are bilineage (reject that they are by sort contamination
% ) in all 3 time points?
n_bilineage = nnz(all(FDRvals_notcontam < 0.05,2));
bilineage = zeros(size(pvals_notcontam,1),max(timepoints));
for t = 1:max(timepoints)
    bilineage(:,t) = all(FDRvals_notcontam(:,timepoints==t) < 0.05,2);
end
is_bilineage = any(bilineage,2);
n_bilineage = nnz(is_bilineage);

%how many clones are monolineage (reject that they are by sort
%contamination in only 1 state in all 3 time points. Same state obviously.)
states = cellfun(@(x) {x(1)}, reads.Properties.VariableNames);
st = {'E','M'};
monolineage = zeros(size(pvals_notcontam,1),2);
for s = 1:2
    FDR = FDRvals_notcontam(:,strcmp(states,st{s})); %FDRs of state
    nFDR = FDRvals_notcontam(:,~strcmp(states,st{s})); %FDRs of other state
    monolineage(:,s) = all(FDR<0.05,2) & ~any(nFDR<0.05,2);
end
if any((monolineage(:,1) & monolineage(:,2)))
    error('something is wrong in monolineage calculation')
end
if any((monolineage(:,1) | monolineage(:,2)) & (is_bilineage))
    error('clones are being marked as mono and bilineage!')
end

n_monolineage = nnz(monolineage(:,1) | monolineage(:,2));

%% test if clones bias change over time

% Let us simply ask, what is the distribution of differences in equilibria
% between time points?

% load in analyzed / normalized clone abundances (post
% contamination-fixing);

rnorm_i = readtable('norm_filtered_read.csv','ReadVariableNames',0);
bc = table2array(readtable('norm_filtered_barcodes.csv','ReadVariableNames',0));
labeli = cellfun(@(x) {x(1:2)},reads_i.Properties.VariableNames);
label = unique(labeli);
rnorm = array2table(zeros(height(rnorm_i),length(label)),'RowNames',bc,'VariableNames',label);
for t = 1:length(label)
    rnorm{:,label{t}} = sum(rnorm_i{:,strcmp(labeli,label{t})},2); %sum postsort replicates
    rnorm{isnan(rnorm{:,label{t}}),label{t}} = 0; %replace NaNs (not found) with 0s
end

%calculate fraction epithelial
timepoints = cellfun(@str2num,(cellfun(@(x) {x(2)}, rnorm.Properties.VariableNames)));
states = cellfun(@(x) {x(1)}, rnorm.Properties.VariableNames);
frac_ep = zeros(height(rnorm),max(timepoints));
for t = 1:max(timepoints)
    E = rnorm{:,timepoints==t & strcmp(states,'E')};
    C = sum(rnorm{:,timepoints==t},2);
    frac_ep(:,t) = E./C;
end

% first calculate the differences by percentage points of percent
% epithelial between 1st and 3rd time points
diffs = 100.*abs(diff(frac_ep(:,[1,3]),1,2));
diffsLessThan15 = nnz(diffs<15)/length(diffs); %fraction of diffs <15%

%% clones mesenchymally biased
mes_bias = nnz(frac_ep(:,1)<.5)/length(frac_ep);


%% calculate shannon entropy of Log2E/M
%calculate log2(E/M)states = cellfun(@(x) {x(1)}, reads.Properties.VariableNames);
states = cellfun(@(x) {x(1)}, rnorm.Properties.VariableNames);
EM = rnorm{:,strcmp(states,'E')}./ rnorm{:,strcmp(states,'M')};
L2EM = (1/3) * sum(log2(EM),2);
%H(X) = - sum(P(xi)*log2(P(xi)));
mes_bias = nnz(L2EM<0)/length(L2EM);
[N,edges] = histcounts(L2EM,min(L2EM):1:max(L2EM));  %bin by 1 and calculate entropy of bins
shannon = N./sum(N).* log2(N./sum(N));
shannon(isnan(shannon)) = 0;
shannon = -sum(shannon);

%% Determine p value of correlation from shuffling clones across time points
L2EM_t = log2(EM);
reps = 1000000;
correls = zeros(reps,size(L2EM_t,2));
n = size(L2EM_t,1);
for i = 1:reps
    fake_eq = L2EM_t(randperm(n),:);
    L = 0;
    for t = 1:size(L2EM_t,2)-1
        for m = 2:size(L2EM_t,2)
            if m <= t
                continue
            end
            L=L+1;
            correls(i,L) = corr(fake_eq(:,m),L2EM_t(:,t));
        end
    end
end

pvals_stability = zeros(size(L2EM_t,2),1);
L = 0;
for t = 1:size(L2EM_t,2)-1
    for m = 2:size(L2EM_t,2)
        if m <= t
            continue
        end
        L=L+1;
        obs = corr(L2EM_t(:,m),L2EM_t(:,t));
        pvals_stability(L) = nnz(correls > obs) / size(correls,1);
    end
end

