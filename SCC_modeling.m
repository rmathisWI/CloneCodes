L2EM = readtable('scclog2EM.csv','ReadVariableNames',0);
L2EM.Properties.VariableNames = {'SCC','Pool'};

reps = 1000000;
correls = zeros(reps,1);
pool = L2EM.Pool;
SCC = L2EM.SCC;
n = length(pool);
for i = 1:reps
    fake_pooleq = pool(randperm(n));
    rho = corrcoef(fake_pooleq,SCC);
    correls(i) = rho(2);
end

[h,edgs] = histcounts(correls,27);
hnrm = h ./ sum(h);
bw = diff(edgs);
cntrs = edgs(2:end) - (bw(1)/2);

res = horzcat(cntrs',hnrm');
nnz(correls > 0.869130933);
writetable(array2table(correls),'random_correlations2.csv');