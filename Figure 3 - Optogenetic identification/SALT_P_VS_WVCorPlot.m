oriList2 = [];

for neuron = 1:numel(neuronIndex.PyrIx_OI)
    if isnan(All{25,1}(neuronIndex.PyrIx_OI(neuron))) | ~ismember(All{25,1}(neuronIndex.PyrIx_OI(neuron)),oriList2)
        oriList2 = [oriList2 neuronIndex.PyrIx_OI(neuron)];
    end
end


n=1;
for cluster2 = 1:numel(oriList2)
        lightIxs = find(All{1,1}{oriList2(cluster2),1}(:,3) < 0.01 & All{1,1}{oriList2(cluster2),1}(:,15)>= parameters_OI.PWCor & All{1,1}{oriList2(cluster2),1}(:,2)>= parameters_OI.minLightInt & All{1,1}{oriList2(cluster2),1}(:,6)>= parameters_OI.ReliabilityTreshold & All{1,1}{oriList2(cluster2),1}(:,4) < parameters_OI.LatencyTreshold & All{1,1}{oriList2(cluster2),1}(:,5) < parameters_OI.JitterTreshold & numel(find(All{24,1}{oriList2(cluster2),1}(:,5) == 1)) > 0);
        [~, stimBlocksDur] = max(All{1,1}{oriList2(cluster2),1}(lightIxs,1));
        [~, stimBlocksInt] = max(All{1,1}{oriList2(cluster2),1}(lightIxs(stimBlocksDur),2));
        try
        latencyTEMP(n) = All{1,1}{oriList2(cluster2),1}(lightIxs(stimBlocksDur(stimBlocksInt)),4);
        JitterTEMP(n) = All{1,1}{oriList2(cluster2),1}(lightIxs(stimBlocksDur(stimBlocksInt)),5);
        ReliabilityTEMP(n) = All{1,1}{oriList2(cluster2),1}(lightIxs(stimBlocksDur(stimBlocksInt)),6);
        WvCorTEMP(n) = All{1,1}{oriList2(cluster2),1}(lightIxs(stimBlocksDur(stimBlocksInt)),15);
        SALTPTEMP(n) = All{1,1}{oriList2(cluster2),1}(lightIxs(stimBlocksDur(stimBlocksInt)),3);
        catch
            [~, stimBlocksDur] = max(All{1,1}{oriList2(cluster2),1}(:,1));
            [~, stimBlocksInt] = max(All{1,1}{oriList2(cluster2),1}(stimBlocksDur,2));
            WvCorTEMP(n) = All{1,1}{oriList2(cluster2),1}(stimBlocksDur(stimBlocksInt),15);
            SALTPTEMP(n) = All{1,1}{oriList2(cluster2),1}(stimBlocksDur(stimBlocksInt),3);
        end
        n=n+1;


end

oriList2 = [];

for neuron = 1:numel(neuronIndex.PyrIx)
    if isnan(All{25,1}(neuronIndex.PyrIx(neuron))) | ~ismember(All{25,1}(neuronIndex.PyrIx(neuron)),oriList2)
        oriList2 = [oriList2 neuronIndex.PyrIx(neuron)];
    end
end


n=1;
for cluster2 = 1:numel(oriList2)
    %lightIxs = find(All{1,1}{cluster2,1}(:,3) < 0.01 & All{1,1}{cluster2,1}(:,15)>= parameters_OI.PWCor & All{1,1}{cluster2,1}(:,2)>= parameters_OI.minLightInt & All{1,1}{cluster2,1}(:,6)>= parameters_OI.ReliabilityTreshold & All{1,1}{cluster2,1}(:,4) < parameters_OI.LatencyTreshold & All{1,1}{cluster2,1}(:,5) < parameters_OI.JitterTreshold & numel(find(All{24,1}{cluster2,1}(:,5) == 1)) > 0);
    [~, stimBlocksDur] = max(All{1,1}{oriList2(cluster2),1}(:,1));
    [~, stimBlocksInt] = max(All{1,1}{oriList2(cluster2),1}(stimBlocksDur,2));
    %latencyTEMP(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),4);
    %JitterTEMP(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),5);
    %ReliabilityTEMP(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),6);
    WvCorTEMP_nonOI(n) = All{1,1}{oriList2(cluster2),1}(stimBlocksDur(stimBlocksInt),15);
    SALTPTEMP_nonOI(n) = All{1,1}{oriList2(cluster2),1}(stimBlocksDur(stimBlocksInt),3);
%     scatter3(All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),4), All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),5), All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),6),60,'d','MarkerFaceColor',[1 0 0])
%     hold on
    n=n+1;
end

figure
hold on
histogram(SALTPTEMP_nonOI,0:0.001:1)
histogram(SALTPTEMP,0:0.001:1)
%scatter(SALTPTEMP_nonOI,WvCorTEMP_nonOI)
%scatter(SALTPTEMP,WvCorTEMP)
xlabel('SALT p-value')
ylabel('WV Cor')
%set(gca,'XScale','log')