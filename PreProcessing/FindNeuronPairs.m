%% Find pairs of neurons that are possibly the same across VAR_ITI and VAR_SD sessions
for neuron = 1:size(neuronIndex.ClustIxPyr,1)
    AnimalNeurons = All{8, 1}{neuronIndex.ClustIxPyr(neuron),1};
    DepthNeurons = All{11, 1}(neuronIndex.ClustIxPyr(neuron),1);
    shankNeurons = All{4, 1}(neuronIndex.ClustIxPyr(neuron),1);
    dateNeurons = All{9, 1}{neuronIndex.ClustIxPyr(neuron),1};
    protNeurons = All{5, 1}{neuronIndex.ClustIxPyr(neuron), 1};
    peakCHNeurons = All{15, 1}(neuronIndex.ClustIxPyr(neuron), 1);
    
    targetNeurons = [];
    for i = 1:size(neuronIndex.ClustIxPyr,1)
        if All{8, 1}{neuronIndex.ClustIxPyr(i),1} == AnimalNeurons & All{11, 1}(neuronIndex.ClustIxPyr(i),1) == DepthNeurons & All{4, 1}(neuronIndex.ClustIxPyr(i),1) == shankNeurons & strcmp(All{9, 1}{neuronIndex.ClustIxPyr(i),1},dateNeurons) == 0 & strcmp(All{5, 1}{neuronIndex.ClustIxPyr(i),1},protNeurons) == 0 & (All{15,1}(neuronIndex.ClustIxPyr(i),1)>=peakCHNeurons-1 & All{15,1}(neuronIndex.ClustIxPyr(i),1)<=peakCHNeurons+1)
            targetNeurons = [targetNeurons neuronIndex.ClustIxPyr(i)];
        end
    end
    
    WVcor = [];
    for partnerNeuron = 1:numel(targetNeurons)
        R = corrcoef(All{13, 1}(neuronIndex.ClustIxPyr(neuron),:),All{13, 1}(targetNeurons(partnerNeuron),:));
        WVcor(partnerNeuron) = R(1,2);
    end
    
    WvCorIx = find(WVcor>0.95);
    ISIcor = [];
    for i = 1:numel(WvCorIx)
        spikeIntervalsTarget = diff(All{21,1}{targetNeurons(WvCorIx(i)),1}(1:All{3, 1}(targetNeurons(WvCorIx(i)),3)))*1000;
        spikeIntervalsNeuron = diff(All{21,1}{neuronIndex.ClustIxPyr(neuron),1}(1:All{3, 1}(neuronIndex.ClustIxPyr(neuron),3)))*1000;
        binSize_hist = 0.5;
        x = [0:binSize_hist:100];
        targetNeuronHist = histcounts(spikeIntervalsTarget,x);
        NeuronHist = histcounts(spikeIntervalsNeuron,x);
        R = corrcoef(targetNeuronHist, NeuronHist);
        ISIcor(i) = R(1,2);
    end
    
    if max(ISIcor)>0.70
        [C I] = max(ISIcor);
        All{25,1}(neuronIndex.ClustIxPyr(neuron),:) = targetNeurons(WvCorIx(I));
    else
        All{25,1}(neuronIndex.ClustIxPyr(neuron),:) = nan;
    end
    
    if isempty(partnerNeuron)
        All{25,1}(neuronIndex.ClustIxPyr(neuron),:) = nan;
    end
    
    clear targetNeurons parterNeuron R
    
end

All{25,1}(All{25,1} == 0) = nan;

% Search for partner neurons that were selected twice. Run correlation over
% PSTH around treshold in that case.
n = [];
for neuron = 1:size(neuronIndex.ClustIxPyr,1)
    if ~isnan(All{25,1}(neuronIndex.ClustIxPyr(neuron))) && numel(find(All{25,1} == All{25,1}(neuronIndex.ClustIxPyr(neuron))))>1
        DoubleIxs = find(All{25,1} == All{25,1}(neuronIndex.ClustIxPyr(neuron)));
        targetNeuron = All{25,1}(DoubleIxs(1));
        protNeurons = All{5, 1}{neuronIndex.ClustIxPyr(neuron), 1};
        if contains(protNeurons,'VAR_SD')
            [~, PSTH_neuron, ~] = PSTH_event_SD(All, neuronIndex, parameters_OI, DoubleIxs(1), 0.2, 'All', 'no');
            [~, PSTH_neuron_two, ~] = PSTH_event_SD(All, neuronIndex, parameters_OI, DoubleIxs(2), 0.2, 'All', 'no');
            [~, PSTH_targetNeuron, ~] = PSTH_event_ITI(All, neuronIndex, parameters_OI, targetNeuron, 0.2, 'All', 'no');
            R_neuron = corrcoef(PSTH_neuron.Cor.wait_start(6:25), PSTH_targetNeuron.Cor.wait_start(6:25));
            R_neuron_two = corrcoef(PSTH_neuron_two.Cor.wait_start(6:25), PSTH_targetNeuron.Cor.wait_start(6:25));
            [~, I] = max([R_neuron(1,2) R_neuron_two(1,2)]);
            DoubleIxs(I) = [];
            All{25,1}(DoubleIxs,1) = nan;
        elseif contains(protNeurons,'VAR_ITI')
            [~, PSTH_neuron, ~] = PSTH_event_ITI(All, neuronIndex, parameters_OI, DoubleIxs(1), 0.2, 'All', 'no');
            [~, PSTH_neuron_two, ~] = PSTH_event_ITI(All, neuronIndex, parameters_OI, DoubleIxs(2), 0.2, 'All', 'no');
            [~, PSTH_targetNeuron, ~] = PSTH_event_SD(All, neuronIndex, parameters_OI, targetNeuron, 0.2, 'All', 'no');
            R_neuron = corrcoef(PSTH_neuron.Cor.wait_start(6:25), PSTH_targetNeuron.Cor.wait_start(6:25));
            R_neuron_two = corrcoef(PSTH_neuron_two.Cor.wait_start(6:25), PSTH_targetNeuron.Cor.wait_start(6:25));
            [~, I] = max([R_neuron(1,2) R_neuron_two(1,2)]);
            DoubleIxs(I) = [];
            All{25,1}(DoubleIxs,1) = nan;
        end
    end
end

    