% Collects all OI spike times and plots them in colormap for Figure 3C

n = 1;
bins_ITI = [];
for neuron = [neuronIndex.PyrIx_VAR_ITI_OI]
try
[~, ~, ~, bins] = PSTH_event_ITI_CorOnlyForFigure(All, neuronIndex, parameters_OI, neuron, 0.2, 'All', 'no','Plot');
bins_ITI(n,:) = bins;
catch
bins_ITI(n,:) = nan;
end
n = n+1;
close all
end

banned = [];
for neuron = 1:numel(neuronIndex.PyrIx_OI)
    if ~isnan(All{25,1}(neuronIndex.PyrIx_OI(neuron)))
        banned = [banned All{25,1}(neuronIndex.PyrIx_OI(neuron))];
    end
end

lists = neuronIndex.PyrIx_OI;
lists(ismember(lists,banned)) = [];

n = 1;
bins_SD = [];
for neuron = [neuronIndex.PyrIx_VAR_SD_OI]
    if ismember(neuron, lists)
        try
        [~, ~, ~, bins] = PSTH_event_SD_CorOnlyForFigure(All, neuronIndex, parameters_OI, neuron, 0.2, 'All', 'no','Plot');
        bins_SD(n,:) = bins;
        catch
        bins_SD(n,:) = nan;
        end
        n = n+1;
        close all
    end
end

%%
[X C] = sort(neuronIndex.PyrIx_VAR_ITI_OI,'ascend');
binSize = 0.001;
figure
colormap('hot')
imagesc(-40:1:40,1:(size([bins_ITI(C,:)],1)),[bins_ITI(C,:)])
axis([-40-(0.5*binSize) 40+(0.5*binSize) 0.5 size([bins_ITI(C,:)],1)])
ylabel('Neuron')
colorbar

[X C] = sort(neuronIndex.PyrIx_VAR_SD_OI,'ascend');
binSize = 0.001;
figure
colormap('hot')
imagesc(-40:1:40,1:(size([bins_SD(C,:)],1)),[bins_SD(C,:)])
axis([-40-(0.5*binSize) 40+(0.5*binSize) 0.5 size([bins_SD(C,:)],1)])
ylabel('Neuron')
colorbar
%%

for i = 1:size(bins_ITI,1)
[C I] = max(bins_ITI(i,:));
bins_ITI(i,:) = bins_ITI(i,:)/C;
end
for i = 1:size(bins_SD,1)
[C I] = max(bins_SD(i,:));
bins_SD(i,:) = bins_SD(i,:)/C;
end


binSize = 0.001;
figure
colormap('hot')
imagesc(-40:1:40,1:(size([bins_SD; bins_ITI],1)),[bins_SD; bins_ITI])
axis([-40-(0.5*binSize) 40+(0.5*binSize) 0.5 size([bins_SD; bins_ITI],1)])
ylabel('Neuron')
colorbar

oriList = neuronIndex.PyrIx_OI;
oriListTEMP = neuronIndex.PyrIx_OI;
oriList2 = [];
banned = [];
for neuron = 1:numel(neuronIndex.PyrIx_OI)
    if isnan(All{25,1}(neuronIndex.PyrIx_OI(neuron))) | ~ismember(All{25,1}(neuronIndex.PyrIx_OI(neuron)),oriList2)
        oriList2 = [oriList2 neuronIndex.PyrIx_OI(neuron)];
        %banned = [banned All{25,1}(neuronIndex.PyrIx_OI(neuron))];
    end
end

%lists = neuronIndex.PyrIx_OI;
%lists(ismember(lists,banned)) = [];

n = 1;
m = 1;
bins2 = [];
bins_SD = [];
bins_ITI = [];
for neuron = 1:numel(oriList2)
        if ismember(oriList2(neuron),neuronIndex.PyrIx_VAR_SD_OI)
            try
            [~, ~, ~, bins] = PSTH_event_SD_CorOnlyForFigure(All, neuronIndex, parameters_OI, oriList2(neuron), 0.2, 'All', 'no','Plot');
            bins2(n,:) = bins;
            catch
            bins2(n,:) = nan;
            end
            n = n+1;
            close all
        elseif ismember(oriList2(neuron),neuronIndex.PyrIx_VAR_ITI_OI)
            try
            [~, ~, ~, bins] = PSTH_event_ITI_CorOnlyForFigure(All, neuronIndex, parameters_OI, oriList2(neuron), 0.2, 'All', 'no','Plot');
            bins2(n,:) = bins;
            catch
            bins2(n,:) = nan;
            end
            n = n+1;
            close all
        end
end

for i = 1:size(bins2,1)
[C I] = max(bins2(i,:));
bins2(i,:) = bins2(i,:)/C;
end

binSize = 0.001;
figure
colormap('hot')
imagesc(-40:1:40,1:(size(bins2,1)),bins2)
axis([-40-(0.5*binSize) 40+(0.5*binSize) 0.5 size(bins2,1)])
ylabel('Neuron')
colorbar
    
