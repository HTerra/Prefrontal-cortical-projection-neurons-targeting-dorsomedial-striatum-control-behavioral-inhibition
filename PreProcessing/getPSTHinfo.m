  

m=1;
ev = {'All','ITI125','ITI50'};
for neuron = 1:numel(neuronIndex.PyrIx_VAR_ITI_OI)
for ev2 = 1
    try
[ITI_OI.(ev{ev2}).psthBinsValue{1,m}, ITI_OI.(ev{ev2}).PSTH{1,m}, ITI_OI.(ev{ev2}).TrialDist{1,m}] = PSTH_event_ITI(All, neuronIndex, parameters_OI, neuronIndex.PyrIx_VAR_ITI_OI(neuron), 0.2, 'All','no');
    catch
    end
    end
m = m+1;
end
save('ITI_OI3.mat','ITI_OI','-v7.3')
clear ITI_OI m

n = 1;
ev = {'All'};
for neuron = 1:numel(neuronIndex.PyrIx_VAR_SD_OI)
for ev2 = 1
    try
[SD_OI.(ev{ev2}).psthBinsValue{1,n}, SD_OI.(ev{ev2}).PSTH{1,n}, SD_OI.(ev{ev2}).TrialDist{1,n}] = PSTH_event_SD(All, neuronIndex, parameters_OI, neuronIndex.PyrIx_VAR_SD_OI(neuron), 0.2, ev{ev2},'no');
    catch
    end
    end
n = n+1;
end
save('SD_OI.mat','SD_OI','-v7.3')
clear SD_OI n

m=1;
ev = {'All'};
for neuron = 1:numel(neuronIndex.PyrIxVAR_ITI)
for ev2 = 1
try
    [ITI_All.(ev{ev2}).psthBinsValue{1,m}, ITI_All.(ev{ev2}).PSTH{1,m}, ITI_All.(ev{ev2}).TrialDist{1,m}] = PSTH_event_ITI(All, neuronIndex, parameters_OI, neuronIndex.PyrIxVAR_ITI(neuron), 0.2, ev{ev2},'no');
catch
end
end
m = m+1;
end
save('ITI_All.mat','ITI_All','-v7.3')
clear ITI_All m

n = 1;
ev = {'All'};
for neuron = 1:numel(neuronIndex.PyrIxVAR_SD)
for ev2 = 1
    try
[SD_All.(ev{ev2}).psthBinsValue{1,n}, SD_All.(ev{ev2}).PSTH{1,n}, SD_All.(ev{ev2}).TrialDist{1,n}] = PSTH_event_SD(All, neuronIndex, parameters_OI, neuronIndex.PyrIxVAR_SD(neuron), 0.2, ev{ev2},'no');
    catch
    end
    end
n = n+1;
end
save('SD_All.mat','SD_All','-v7.3')
clear SD_All n