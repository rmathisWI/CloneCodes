function DensityPlot(data1,data2)
[values, centers] = hist3([data1(~isinf(data2) & ~isinf(data1)) data2(~isinf(data2) & ~isinf(data1))],[64 64]);
im = imagesc(centers{:}, values');
alphamask = values'>0;
im.AlphaData = alphamask;

end

