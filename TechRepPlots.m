reads = readtable('norm_filtered_read.csv','ReadVariableNames',0);
reads.Properties.VariableNames = {'E1_1','E1_2','M1_1','M1_2','E2_1','E2_2','M2_1','M2_2','E3_1','E3_2','M3_1','M3_2'};
timepoints = cellfun(@str2num,(cellfun(@(x) {x(2)}, reads.Properties.VariableNames)));
states = cellfun(@(x) {x(1)}, reads.Properties.VariableNames);

figure()
fig = gcf;
fig.Units='inches'; 
fig.InvertHardcopy='off';
%fig.Position(2:4) = [4 3.8542    6.2187]; 
fig.Position(2:4) = [4 6.2187   3.8542]; 
fig.PaperPositionMode='Auto';
fig.Renderer='painters';
n=0;
rhos = zeros(max(timepoints),2);
labels = {'Day 0', 'Day 7', 'Day 14'};
for s = 1:2
    for t = 1:max(timepoints);
        tp = timepoints==t;
        n=n+1;
        switch s
            case 1
                st = 'E';
                ttl = 'Epi. ';
            case 2
                st = 'M';
                ttl = 'Mes. ';
        end
        sts = strcmp(states,st);
        reads_t = reads{:,tp & sts};
        lbls = reads.Properties.VariableNames(tp & sts);
        
        %subplot(max(timepoints),2,n)
        subplot(2,max(timepoints),n)
        
        DensityPlot(log10(reads_t(:,1)),log10(reads_t(:,2)));
        ax = gca;
        ax.YDir = 'normal';
        ax.XLim = [-8 -1];
        ax.YLim = ax.XLim;
        ax.Color = 'none';
        ax.Box = 'off';
        ax.TickDir = 'out';
        ax.FontSize = 10;
        ax.LineWidth=1;
        ax.TitleFontSizeMultiplier = 1;
        ax.LabelFontSizeMultiplier = 1;
        ax.FontName = 'Arial';
        title(horzcat(ttl,labels{t}));
        ylabel({'Log_1_0 Fraction Repl. 2)'});
        xlabel({'Log_1_0 Fraction Repl. 1)'});
        r = corrcoef(log10(reads_t)); %pearson correlation coefficient
        %r = round(r(2),2,'Significant');
        r = r(2);
        text(0.05,0.95,horzcat('r = ',num2str(r,3)),'Units','normalized','FontSize',10,'FontName','Arial');
        rhos(t,s) = r;
    end
end
fig.Color = 'none';
print('TechReps.svg','-dsvg','-painters');

        