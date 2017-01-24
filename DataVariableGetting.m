%get data
reads_i = readtable('norm_filtered_read - Copy.csv','ReadRowNames',1,'ReadVariableNames',1,'HeaderLines',1); %load in reads of clones
reads_i.Properties.VariableNames = {'E1_1','E1_2','M1_1','M1_2','E2_1','E2_2','M2_1','M2_2','E3_1','E3_2','M3_1','M3_2'};

%sum postsort replicates
labeli = cellfun(@(x) {x(1:2)},reads_i.Properties.VariableNames);
label = unique(labeli);
reads = array2table(zeros(height(reads_i),length(label)),'RowNames',reads_i.Properties.RowNames,'VariableNames',label);
for t = 1:length(label)
    reads{:,label{t}} = sum(reads_i{:,strcmp(labeli,label{t})},2); %sum postsort replicates
    reads{isnan(reads{:,label{t}}),label{t}} = 0; %replace NaNs (not found) with 0s
end
reads2 = reads;
reads2{:,:} = reads{:,:}./2; %make fraction of total population


%% first: distribution of clones Log2EM 
states = cellfun(@(x) {x(1)}, reads.Properties.VariableNames);
EM = reads{:,strcmp(states,'E')}./ reads{:,strcmp(states,'M')};
L2EM = (1/3) * sum(log2(EM),2);


%% Second: growth rate of epithelial and mesenchymal clones

% estimate parameters of population

td = 3; %estimated population doublings in 1 week
% N = N0 e^(k(t1 - t0) ; k = ln(N/N0)/(t1-t0)
N = 1*2^td; %N if N0 =1;
k = log(N)/7; %where 7 is days

% estimate population sizes at day 0, 7, 14;
N0 = 1;
t1 = 7;
t2 = 14;
N1 = N0*exp(k*t1); % N = N0 e^(k(t1 - t0)
N2 = N0*exp(k*t2); % N = N0 e^(k(t1 - t0)
Ns = vertcat(N0,N1,N2);

% now estimate cell numbers for clones
cellNs = zeros(size(reads2)); %instantiate matrix for cell numbers for each clone
timepoints = cellfun(@str2num,(cellfun(@(x) {x(2)}, reads2.Properties.VariableNames)));
for t = 1:max(timepoints) %for each timepoint
    cellNs(:,timepoints==t) = reads2{:,timepoints==t} .* Ns(t); %multiply fraction by number of cells
end

% compute clone growth rates
%k = ln(N/N0)/(t1-t0)
ks = zeros(height(reads2),max(timepoints)-1); %rate constant matrix
for t = 2 : max(timepoints) %for each timepoint from first
    %ks(:,t-1) = log( sum(cellNs(:,timepoints==t),2) ./ sum(cellNs(:,timepoints == t-1),2)) ./ 7;
    ks(:,t-1) = log( sum(cellNs(:,timepoints==t),2) ./ sum(cellNs(:,timepoints == 1),2)) ./ (7*(t-1));
end
k_mean = mean(ks,2);

%% record in csv
observeddata = table(k_mean,L2EM,sum(reads2{:,timepoints==1}.*2.9E7,2),'VariableNames',{'k','L2EM','cellNs'});
writetable(observeddata,'observeddata.csv','WriteVariableNames',1);



figure()
fig = gcf;
fig.InvertHardcopy='off';
fig.Units='inches'; 
fig.Position(3:4) = [3    2.5];
fig.PaperPositionMode='Auto';
fig.Renderer='painters';
fig.Color = 'none';
DensityPlot(ks(:,1),ks(:,2));
ax=gca;
ax.YDir = 'Normal';
ax.Box = 'off';
ax.TickDir = 'out';
ax.Color = 'none';
ax.FontName = 'Arial';
ax.FontSize = 10;
ax.LineWidth=1;
axs = [ax.XLim;ax.YLim];
ax.XLim(1) = min(axs(:,1));
ax.XLim(2) = max(axs(:,2));
ax.YLim = ax.XLim;
ax.YTick = ax.XTick;
ax.TitleFontSizeMultiplier=1;
ax.LabelFontSizeMultiplier = 1; 
X = ax.XLim;
line(X,X,'LineWidth',1.5,'Color','k','LineStyle','--');
xlabel({'k /day','Day 0 to Day 7'});
ylabel({'k /day','Day 0 to Day 14'});
text(0.1,0.9,horzcat('rho = ',num2str(corr(ks(:,1),ks(:,2)))),'units','normalized','FontName','Arial','FontSize',10);
print('growthrates.svg','-dsvg');

figure()
fig = gcf;
fig.InvertHardcopy='off';
fig.Units='inches'; 
fig.Position(3:4) = [3    2.5];
fig.PaperPositionMode='Auto';
fig.Renderer='painters';
fig.Color = 'none';
h=histogram(k_mean);
h.FaceColor = 'none';
h.LineWidth = 1;
ax=gca;
ax.YDir = 'Normal';
ax.Box = 'off';
ax.TickDir = 'out';
ax.Color = 'none';
ax.FontName = 'Arial';
ax.FontSize = 10;
ax.LineWidth=1;
ax.TitleFontSizeMultiplier=1;
ax.LabelFontSizeMultiplier = 1; 
xlabel('Growth Rate (/ day)');
ylabel('# of clones in bin')
print('growthrate hist.svg','-dsvg');



figure()
fig = gcf;
fig.InvertHardcopy='off';
fig.Units='inches'; 
fig.Position(3:4) = [3    2.5];
fig.PaperPositionMode='Auto';
fig.Renderer='painters';
fig.Color = 'none';
scatter(L2EM,k_mean,5,'filled');
ax=gca;
ax.YDir = 'Normal';
ax.Box = 'off';
ax.TickDir = 'out';
ax.Color = 'none';
ax.FontName = 'Arial';
ax.FontSize = 10;
ax.LineWidth=1;
ax.TitleFontSizeMultiplier=1;
ax.LabelFontSizeMultiplier = 1; 
xlabel('log2[E/M]')
ylabel('mean growth rate (/day)');
text(0.1,0.9,horzcat('rho = ',num2str(corr(k_mean,L2EM))),'units','normalized','FontSize',10,'FontName','Arial');
print('growthrate vs EM.svg','-dsvg');

figure()
scatter(k_mean,log10(cellNs(:,1)),5,'filled');
xlabel('mean growth rate')
ylabel('log10 cell numbers day0');
text(0.1,0.9,horzcat('rho = ',num2str(corr(k_mean,cellNs(:,1)))),'units','normalized');

%% are more-epithelial clones larger?
N0 = 2.9E7; %cells at first time point
cellN0 = reads2{:,timepoints==1} .* N0;
figure()
fig = gcf;
fig.InvertHardcopy='off';
fig.Units='inches'; 
fig.Position(3:4) = [3    2.5];
fig.PaperPositionMode='Auto';
fig.Renderer='painters';
fig.Color = 'none';
hold on
[y,x] = ecdf(cellN0(L2EM<0));
plot(x,y,'LineWidth',3);
[y,x] = ecdf(cellN0(L2EM>0));
plot(x,y,'LineWidth',3);
text(0.1,0.8,horzcat('p = ',num2str(ranksum(cellN0(L2EM<0),cellN0(L2EM>0)))),'units','normalized','FontName','Arial','FontSize',10);
ax = gca;
ax.YDir = 'Normal';
ax.Box = 'off';
ax.TickDir = 'out';
ax.Color = 'none';
ax.FontName = 'Arial';
ax.FontSize = 10;
ax.LineWidth=1;
ax.TitleFontSizeMultiplier=1;
ax.LabelFontSizeMultiplier = 1; 
ax.XScale = 'log';
ax.XLim =[1 1E6];
ax.XTick = [1 100  10000 1000000];
L = legend({'more mes. clones','more epi. clones'},'location','NW','FontSize',10,'FontName','Arial');
L.Box = 'off';
xlabel('Cell Numbers at Day 0')
ylabel('Cumulative Density');
print('clonesize vs EM.svg','-dsvg');

%% estimate growth rate distributions of more E and more M clones
GRdistE = fitdist(k_mean(L2EM>0),'normal');
GRdistM = fitdist(k_mean(L2EM<0),'normal');
