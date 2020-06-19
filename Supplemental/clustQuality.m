function clustQuality(All, ClustIx, OIindex)

% ID plot
temp = figure;
y3ID = histogram(All{3,1}(ClustIx,1),'Normalization','cdf');
hold on
y4ID = histogram(All{3,1}(OIindex,1),'Normalization','cdf');
y3ID.BinWidth = 20;
y4ID.BinWidth = 20;
hold off

figure
hold on
y1ID = histogram(All{3,1}(ClustIx,1),'FaceColor','b', 'LineWidth', 2);
y2ID = histogram(All{3,1}(OIindex,1),'FaceColor','r', 'LineWidth', 2);
y1ID.BinWidth = 20;
y2ID.BinWidth = 20;

cumPlotID = [ 0 y3ID.Values]*max(y1ID.Values);
cumPlotOIID = [0 y4ID.Values 1]*max(y1ID.Values);
cumPlotID(end+1:length(cumPlotOIID)) = max(cumPlotID);

Limit = max([y4ID.BinLimits y3ID.BinLimits]);
xID = linspace(40,Limit,Limit/y1ID.BinWidth-1);
stairs(xID,[cumPlotID(2:end) cumPlotID(end)],'b','LineWidth',2)
stairs(xID,[cumPlotOIID(2:end) cumPlotOIID(end) cumPlotOIID(end) cumPlotOIID(end) cumPlotOIID(end) cumPlotOIID(end) cumPlotOIID(end) cumPlotOIID(end)],'r','LineWidth',2)


title('Cluster Isolation Distance %s')
xlabel('Cluster Isolation Distance')
ylabel('No. of neurons')

close(temp)



% LR plot
temp = figure;
y3LR = histogram(All{3,1}(ClustIx,2),'Normalization','cdf');
hold on
y4LR = histogram(All{3,1}(OIindex,2),'Normalization','cdf');
y3LR.BinWidth = 0.1;
y4LR.BinWidth = 0.1;
hold off

figure
hold on
y1LR = histogram(All{3,1}(ClustIx,2),'FaceColor','b', 'LineWidth', 2);
y2LR = histogram(All{3,1}(OIindex,2),'FaceColor','r', 'LineWidth', 2);
y1LR.BinWidth = 0.1;
y2LR.BinWidth = 0.1;

cumPlotLR = y3LR.Values*max(y1LR.Values);
cumPlotOILR = y4LR.Values*max(y1LR.Values);
cumPlotLR(end+1:length(cumPlotOILR)) = max(cumPlotLR);

Limit = max([y4LR.BinLimits y3LR.BinLimits]);
xLR = linspace(0,Limit,Limit/y1LR.BinWidth);
stairs(xLR,[cumPlotLR(2:end) cumPlotLR(end)],'b','LineWidth',2)
stairs(xLR,[cumPlotOILR(2:end) cumPlotOILR(end)],'r','LineWidth',2)

title('Cluster L-ratio')
xlabel('Cluster L-ratio')
ylabel('No. of neurons')

close(temp)

% ISI plot
temp = figure;
y3ISI = histogram(All{3,1}(ClustIx,7),'Normalization','cdf');
hold on
y4ISI = histogram(All{3,1}(OIindex,7),'Normalization','cdf');
y3ISI.BinWidth = 0.01;
y4ISI.BinWidth = 0.01;
hold off

figure
hold on
y1ISI = histogram(All{3,1}(ClustIx,7),'FaceColor','b', 'LineWidth', 2);
y2ISI = histogram(All{3,1}(OIindex,7),'FaceColor','r', 'LineWidth', 2);
y1ISI.BinWidth = 0.01;
y2ISI.BinWidth = 0.01;

cumPlotISI = y3ISI.Values*max(y1ISI.Values);
cumPlotOIISI = y4ISI.Values*max(y1ISI.Values);
cumPlotISI(end+1:length(cumPlotOIISI)) = max(cumPlotISI);

Limit = max([y4ISI.BinLimits y3ISI.BinLimits]);
xISI = linspace(0,Limit,Limit/y1ISI.BinWidth);
stairs(xISI,[cumPlotISI(2:end) cumPlotISI(end)],'b','LineWidth',2)
stairs(xISI,[cumPlotOIISI(2:end) cumPlotOIISI(end)],'r','LineWidth',2)

title('Cluster ISI < 1.5ms')
xlabel('Cluster ISI < 1.5ms (fraction)')
ylabel('No. of neurons')

close(temp)