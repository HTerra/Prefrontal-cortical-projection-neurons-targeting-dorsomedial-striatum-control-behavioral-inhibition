% Function plots locations of different neurons categories per different
% session type and neuron group (tagged and non-tagged)

%% Get non-oi neuron Ixs from session that have OI neurons in it

VAR_ITIixOIMatched = [];
for neuron = 1:numel(neuronIndex.PyrIxVAR_ITI)
    for i = neuronIndex.PyrIx_VAR_ITI_OI
        if strcmp(All{9, 1}{neuronIndex.PyrIxVAR_ITI(neuron),1},All{9, 1}{i,1}) & numel(find(All{1, 1}{neuronIndex.PyrIxVAR_ITI(neuron), 1}(:,3)<0.01))==0
            VAR_ITIixOIMatched = [VAR_ITIixOIMatched neuron];
            break
        else
        end
    end
end

VAR_SDixOIMatched = [];
for neuron = 1:numel(neuronIndex.PyrIxVAR_SD)
    for i = neuronIndex.PyrIx_VAR_SD_OI
        if strcmp(All{9, 1}{neuronIndex.PyrIxVAR_SD(neuron),1},All{9, 1}{i,1}) & numel(find(All{1, 1}{neuronIndex.PyrIxVAR_SD(neuron), 1}(:,3)<0.01))==0
            VAR_SDixOIMatched = [VAR_SDixOIMatched neuron];
            break
        else
        end
    end
end

%% Get depth profiles of non-OI and OI neurons over different sessions.

DepthFromPia_SD = [];
DepthFromSurface_SD = [];

actIx = find(SD_All_clust.TH(VAR_SDixOIMatched) == 1);
otherIx = find(SD_All_clust.TH(VAR_SDixOIMatched) == 0);
silIx = find(SD_All_clust.TH(VAR_SDixOIMatched) == -1);

for neuron = 1:size(neuronIndex.PyrIxVAR_SD(VAR_SDixOIMatched),2)
    DepthFromPia_SD = [DepthFromPia_SD All{4, 1}(neuronIndex.PyrIxVAR_SD(VAR_SDixOIMatched(neuron)),2)];
    DepthFromSurface_SD = [DepthFromSurface_SD All{11, 1}(neuronIndex.PyrIxVAR_SD(VAR_SDixOIMatched(neuron)),1)];
end

DepthFromPia_ITI = [];
DepthFromSurface_ITI = [];

for neuron = 1:size(neuronIndex.PyrIxVAR_ITI(VAR_ITIixOIMatched),2)
    DepthFromPia_ITI = [DepthFromPia_ITI All{4, 1}(neuronIndex.PyrIxVAR_ITI(VAR_ITIixOIMatched(neuron)),2)];
    DepthFromSurface_ITI = [DepthFromSurface_ITI All{11, 1}(neuronIndex.PyrIxVAR_ITI(VAR_ITIixOIMatched(neuron)),1)];
end

DepthFromPia_ITIOI = [];
DepthFromSurface_ITIOI = [];

actIx = find(ITI_OI_clust.TH == 1);
otherIx = find(ITI_OI_clust.TH == 0);
silIx = find(ITI_OI_clust.TH == -1);
for neuron = 1:size(neuronIndex.PyrIx_VAR_ITI_OI,2)
    DepthFromPia_ITIOI = [DepthFromPia_ITIOI All{4, 1}(neuronIndex.PyrIx_VAR_ITI_OI(neuron),2)];
    DepthFromSurface_ITIOI = [DepthFromSurface_ITIOI All{11, 1}(neuronIndex.PyrIx_VAR_ITI_OI(neuron),1)];
end


%% Plot

figure
subplot(2,2,1)
hold on
% scatter(DepthFromPia_ITIOI(actIx)+15,DepthFromSurface_ITIOI(actIx),10,'k','filled')
% scatter(DepthFromPia_ITIOI(otherIx),DepthFromSurface_ITIOI(otherIx),10,'k','filled')
% scatter(DepthFromPia_ITIOI(silIx)-15,DepthFromSurface_ITIOI(silIx),10,'k','filled')
scatter(DepthFromPia_ITIOI,DepthFromSurface_ITIOI,10,'k','filled')
xlim([0 1400])
ylim([1200 3400])
xlabel('Depth from pia')
ylabel('Depth from surface')
title('var wait time OI')
axis ij
subplot(2,2,2)
hold on
%AllNeuronsBinsDepth = histcounts(DepthFromSurface_ITIOI,[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400])+histcounts(DepthFromSurface_ITI,[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]);
plot([[1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]],histcounts(DepthFromSurface_ITIOI(actIx),[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]),'r');%/(numel(DepthFromPia_ITI)+numel(DepthFromPia_ITIOI)),'r')
plot([[1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]],histcounts(DepthFromSurface_ITIOI(otherIx),[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]),'k');%/(numel(DepthFromPia_ITI)+numel(DepthFromPia_ITIOI)),'k')
plot([[1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]],histcounts(DepthFromSurface_ITIOI(silIx),[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]),'b');%/(numel(DepthFromPia_ITI)+numel(DepthFromPia_ITIOI)),'b')
xlim([1200 3400])
subplot(2,2,3)
hold on
%AllNeuronsBinsPia = histcounts(DepthFromPia_ITIOI,[0 200 400 600 800 1000 1200 1400])+histcounts(DepthFromPia_ITI,[0 200 400 600 800 1000 1200 1400]);
plot([200 400 600 800 1000 1200 1400],histcounts(DepthFromPia_ITIOI(actIx),[0 200 400 600 800 1000 1200 1400]),'r');%/(numel(DepthFromPia_ITI)+numel(DepthFromPia_ITIOI)),'r')
plot([200 400 600 800 1000 1200 1400],histcounts(DepthFromPia_ITIOI(otherIx),[0 200 400 600 800 1000 1200 1400]),'k');%/(numel(DepthFromPia_ITI)+numel(DepthFromPia_ITIOI)),'k')
plot([200 400 600 800 1000 1200 1400],histcounts(DepthFromPia_ITIOI(silIx),[0 200 400 600 800 1000 1200 1400]),'b');%/(numel(DepthFromPia_ITI)+numel(DepthFromPia_ITIOI)),'b')
xlim([0 1400])

DepthFromPia_SDOI = [];
DepthFromSurface_SDOI = [];

actIx = find(SD_OI_clust.TH == 1);
otherIx = find(SD_OI_clust.TH == 0);
silIx = find(SD_OI_clust.TH == -1);

for neuron = 1:size(neuronIndex.PyrIx_VAR_SD_OI,2)
    DepthFromPia_SDOI = [DepthFromPia_SDOI All{4, 1}(neuronIndex.PyrIx_VAR_SD_OI(neuron),2)];
    DepthFromSurface_SDOI = [DepthFromSurface_SDOI All{11, 1}(neuronIndex.PyrIx_VAR_SD_OI(neuron),1)];
end

figure
subplot(2,2,1)
hold on
% scatter(DepthFromPia_SDOI(actIx)+15,DepthFromSurface_SDOI(actIx),10,'k','filled')
% scatter(DepthFromPia_SDOI(otherIx),DepthFromSurface_SDOI(otherIx),10,'k','filled')
% scatter(DepthFromPia_SDOI(silIx)-15,DepthFromSurface_SDOI(silIx),10,'k','filled')
scatter(DepthFromPia_SDOI,DepthFromSurface_SDOI,10,'k','filled')
xlim([0 1400])
ylim([1200 3400])
xlabel('Depth from pia')
ylabel('Depth from surface')
title('var stim duration OI')
axis ij
subplot(2,2,2)
hold on
%AllNeuronsBinsDepth = histcounts(DepthFromSurface_SDOI,[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400])+histcounts(DepthFromSurface_SD,[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]);
plot([[1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]],histcounts(DepthFromSurface_SDOI(actIx),[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]),'r');%/(numel(DepthFromPia_SD)+numel(DepthFromPia_SDOI)),'r')
plot([[1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]],histcounts(DepthFromSurface_SDOI(otherIx),[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]),'k');%/(numel(DepthFromPia_SD)+numel(DepthFromPia_SDOI)),'k')
plot([[1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]],histcounts(DepthFromSurface_SDOI(silIx),[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]),'b');%/(numel(DepthFromPia_SD)+numel(DepthFromPia_SDOI)),'b')
xlim([1200 3400])
subplot(2,2,3)
hold on
%AllNeuronsBinsPia = histcounts(DepthFromPia_ITIOI,[0 200 400 600 800 1000 1200 1400])+histcounts(DepthFromPia_ITI,[0 200 400 600 800 1000 1200 1400]);
plot([200 400 600 800 1000 1200 1400],histcounts(DepthFromPia_SDOI(actIx),[0 200 400 600 800 1000 1200 1400]),'r');%/(numel(DepthFromPia_SD)+numel(DepthFromPia_SDOI)),'r')
plot([200 400 600 800 1000 1200 1400],histcounts(DepthFromPia_SDOI(otherIx),[0 200 400 600 800 1000 1200 1400]),'k');%/(numel(DepthFromPia_SD)+numel(DepthFromPia_SDOI)),'k')
plot([200 400 600 800 1000 1200 1400],histcounts(DepthFromPia_SDOI(silIx),[0 200 400 600 800 1000 1200 1400]),'b');%/(numel(DepthFromPia_SD)+numel(DepthFromPia_SDOI)),'b')
xlim([0 1400])


actIx = find(ITI_All_clust.TH(VAR_ITIixOIMatched) == 1);
otherIx = find(ITI_All_clust.TH(VAR_ITIixOIMatched) == 0);
silIx = find(ITI_All_clust.TH(VAR_ITIixOIMatched) == -1);

figure
subplot(2,2,1)
hold on
% scatter(DepthFromPia_ITI(actIx)+15,DepthFromSurface_ITI(actIx),10,'k','filled')
% scatter(DepthFromPia_ITI(otherIx),DepthFromSurface_ITI(otherIx),10,'k','filled')
% scatter(DepthFromPia_ITI(silIx)-15,DepthFromSurface_ITI(silIx),10,'k','filled')
scatter(DepthFromPia_ITI,DepthFromSurface_ITI,10,'k','filled')
xlim([0 1400])
ylim([1200 3400])
xlabel('Depth from pia')
ylabel('Depth from surface')
title('var wait time non-OI')
axis ij
subplot(2,2,2)
hold on
plot([[1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]],histcounts(DepthFromSurface_ITI(actIx),[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]),'r');%/(numel(DepthFromPia_ITI)+numel(DepthFromPia_ITIOI)),'r')
plot([[1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]],histcounts(DepthFromSurface_ITI(otherIx),[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]),'k');%/(numel(DepthFromPia_ITI)+numel(DepthFromPia_ITIOI)),'k')
plot([[1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]],histcounts(DepthFromSurface_ITI(silIx),[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]),'b');%/(numel(DepthFromPia_ITI)+numel(DepthFromPia_ITIOI)),'b')
xlim([1200 3400])
subplot(2,2,3)
hold on
plot([200 400 600 800 1000 1200 1400],histcounts(DepthFromPia_ITI(actIx),[0 200 400 600 800 1000 1200 1400]),'r');%/(numel(DepthFromPia_ITI)+numel(DepthFromPia_ITIOI)),'r')
plot([200 400 600 800 1000 1200 1400],histcounts(DepthFromPia_ITI(otherIx),[0 200 400 600 800 1000 1200 1400]),'k');%/(numel(DepthFromPia_ITI)+numel(DepthFromPia_ITIOI)),'k')
plot([200 400 600 800 1000 1200 1400],histcounts(DepthFromPia_ITI(silIx),[0 200 400 600 800 1000 1200 1400]),'b');%/(numel(DepthFromPia_ITI)+numel(DepthFromPia_ITIOI)),'b')
xlim([0 1400])

actIx = find(SD_All_clust.TH(VAR_SDixOIMatched) == 1);
otherIx = find(SD_All_clust.TH(VAR_SDixOIMatched) == 0);
silIx = find(SD_All_clust.TH(VAR_SDixOIMatched) == -1);

figure
subplot(2,2,1)
hold on
%scatter(DepthFromPia_SD(actIx)+15,DepthFromSurface_SD(actIx),10,'k','filled')
%scatter(DepthFromPia_SD(otherIx),DepthFromSurface_SD(otherIx),10,'k','filled')
%scatter(DepthFromPia_SD(silIx)-15,DepthFromSurface_SD(silIx),10,'k','filled')
scatter(DepthFromPia_SD,DepthFromSurface_SD,10,'k','filled')
xlim([0 1400])
ylim([1200 3400])
xlabel('Depth from pia')
ylabel('Depth from surface')
title('var stim duration non-OI')
axis ij
subplot(2,2,2)
hold on
plot([[1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]],histcounts(DepthFromSurface_SD(actIx),[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]),'r');%/(numel(DepthFromPia_SD)+numel(DepthFromPia_SDOI)),'r')
plot([[1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]],histcounts(DepthFromSurface_SD(otherIx),[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]),'k');%/(numel(DepthFromPia_SD)+numel(DepthFromPia_SDOI)),'k')
plot([[1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]],histcounts(DepthFromSurface_SD(silIx),[1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400]),'b');%/(numel(DepthFromPia_SD)+numel(DepthFromPia_SDOI)),'b')
xlim([1200 3400])
subplot(2,2,3)
hold on
plot([200 400 600 800 1000 1200 1400],histcounts(DepthFromPia_SD(actIx),[0 200 400 600 800 1000 1200 1400]),'r');%/(numel(DepthFromPia_SD)+numel(DepthFromPia_SDOI)),'r')
plot([200 400 600 800 1000 1200 1400],histcounts(DepthFromPia_SD(otherIx),[0 200 400 600 800 1000 1200 1400]),'k');%/(numel(DepthFromPia_SD)+numel(DepthFromPia_SDOI)),'k')
plot([200 400 600 800 1000 1200 1400],histcounts(DepthFromPia_SD(silIx),[0 200 400 600 800 1000 1200 1400]),'b');%/(numel(DepthFromPia_SD)+numel(DepthFromPia_SDOI)),'b')
xlim([0 1400])

