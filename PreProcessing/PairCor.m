function PairCor = PairCor(All,neuronIndex)

% Function returns linear (Spearman's) corr coeff between the PSTH 1 sec to
% 3 sec around cue orientation of neuron pairs recorded in VAR-ITI and
% VAR_SD conditions.

Correl = [];
for neuron = 1:size(All{25,1},1)
    if ~isnan(All{25,1}(neuron,1)) && contains(All{5, 1}{neuron, 1},'ITI')
        [~, PSTH, ~] = PSTH_event_ITI(All, neuronIndex, parameters_OI, neuron, 0.2, 'All', 'no');
        [~, PSTH2, ~] = PSTH_event_SD(All, neuronIndex, parameters_OI, All{25,1}(neuron,1), 0.2, 'All', 'no');
        R = corrcoef(PSTH.Cor.wait_start(6:25), PSTH2.Cor.wait_start(6:25));
        Correl = [Correl R(1,2)];
    elseif ~isnan(All{25,1}(neuron,1)) && contains(All{5, 1}{neuron, 1},'SD')
        [~, PSTH, ~] = PSTH_event_SD(All, neuronIndex, parameters_OI, neuron, 0.2, 'All', 'no');
        [~, PSTH2, ~] = PSTH_event_ITI(All, neuronIndex, parameters_OI, All{25,1}(neuron,1), 0.2, 'All', 'no');
        R = corrcoef(PSTH.Cor.wait_start(6:25), PSTH2.Cor.wait_start(6:25));
        Correl = [Correl R(1,2)];
    end
end

figure
boxplot(Correl)
ylabel('Correlation coefficient')

peak1 = [];
peak2 = [];
for neuron = 1:size(All{25,1},1)
    if ~isnan(All{25,1}(neuron,1)) && contains(All{5, 1}{neuron, 1},'ITI')
        [~, PSTH, ~] = PSTH_event_ITI(All, neuronIndex, parameters_OI, neuron, 0.2, 'All', 'no');
        [~, PSTH2, ~] = PSTH_event_SD(All, neuronIndex, parameters_OI, All{25,1}(neuron,1), 0.2, 'All', 'no');
        peak1 = [peak1 max(abs(PSTH.Cor.wait_start(6:25)))];
        peak2 = [peak2 max(abs(PSTH2.Cor.wait_start(6:25)))];
    elseif ~isnan(All{25,1}(neuron,1)) && contains(All{5, 1}{neuron, 1},'SD')
        [~, PSTH, ~] = PSTH_event_SD(All, neuronIndex, parameters_OI, neuron, 0.2, 'All', 'no');
        [~, PSTH2, ~] = PSTH_event_ITI(All, neuronIndex, parameters_OI, All{25,1}(neuron,1), 0.2, 'All', 'no');
        peak1 = [peak1 max(abs(PSTH.Cor.wait_start(6:25)))];
        peak2 = [peak2 max(abs(PSTH2.Cor.wait_start(6:25)))];
    end
end

figure
scatter(peak1,peak2)