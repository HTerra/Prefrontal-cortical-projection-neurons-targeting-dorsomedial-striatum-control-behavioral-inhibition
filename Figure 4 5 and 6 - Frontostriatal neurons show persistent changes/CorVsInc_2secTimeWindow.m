function CorVSinc_2secTimeWindow(All,neuronIndex, SD_All, SD_OI, SDperiod)

% Makes Figure S4D

binSize = 0.2;
SDperiod = 'All';
trialTypo = {'TS','TH','Cue'};
trialTypo2 = {'TS','TH','Cue'};
timeBins = {11:15, 48:57, 26:35};

%% 1. Get CorFr for OI and All per sync point

for i = 1:numel(neuronIndex.PyrIx_VAR_SD_OI)
    PSTH_Cor_OI{i} = [SD_OI.(SDperiod).psthBinsValue{1, i}.Cor.trial_start(:,1:35) SD_OI.(SDperiod).psthBinsValue{1, i}.Cor.wait_start(:,1:25)];
    PSTH_Inc_OI{i} = [SD_OI.(SDperiod).psthBinsValue{1, i}.Inc.trial_start(:,1:35) SD_OI.(SDperiod).psthBinsValue{1, i}.Inc.wait_start(:,1:25)];
end

for i = 1:numel(neuronIndex.PyrIxVAR_SD)
    PSTH_Cor_All{i} = [SD_All.(SDperiod).psthBinsValue{1, i}.Cor.trial_start(:,1:35) SD_All.(SDperiod).psthBinsValue{1, i}.Cor.wait_start(:,1:25)];
    PSTH_Inc_All{i} = [SD_All.(SDperiod).psthBinsValue{1, i}.Inc.trial_start(:,1:35) SD_All.(SDperiod).psthBinsValue{1, i}.Inc.wait_start(:,1:25)];
end

%% 1. Get CorFr for OI and All per sync point, calculate significance

for t = 1:size(timeBins,2)
    for i = 1:numel(neuronIndex.PyrIx_VAR_SD_OI)
            p = ranksum(mean(PSTH_Cor_OI{1, i}(:,timeBins{t}),2)/binSize,mean(PSTH_Inc_OI{1, i}(:,timeBins{t})/binSize,2));
            Sign.SD_OI.(trialTypo2{t})(i) = p;
            CorFR.OI.(trialTypo2{t})(i) = mean(mean(PSTH_Cor_OI{1, i}(:,timeBins{t})/binSize,2)-mean(PSTH_Cor_OI{1, i}(:,1:10)/binSize,2));
            IncFR.OI.(trialTypo2{t})(i) = mean(mean(PSTH_Inc_OI{1, i}(:,timeBins{t})/binSize,2)-mean(PSTH_Inc_OI{1, i}(:,1:10)/binSize,2));
            CorFR_AUC.OI.(trialTypo2{t})(i,:) = sum(mean(PSTH_Cor_OI{1, i}(:,timeBins{t})/binSize,1)-mean(mean(PSTH_Cor_OI{1, i}(:,1:10)/binSize,1)));
            IncFR_AUC.OI.(trialTypo2{t})(i,:) = sum(mean(PSTH_Inc_OI{1, i}(:,timeBins{t})/binSize,1)-mean(mean(PSTH_Inc_OI{1, i}(:,1:10)/binSize,1)));
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

for t = 1:size(timeBins,2)
    for i = 1:numel(neuronIndex.PyrIxVAR_SD(VAR_SDixOIMatched))
            p = ranksum(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t}),2)/binSize,mean(PSTH_Inc_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t}),2)/binSize);
            Sign.SD_All.(trialTypo2{t})(i) = p;
            CorFR.All.(trialTypo2{t})(i) = mean(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t})/binSize,2)-mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,1:10)/binSize,2));
            IncFR.All.(trialTypo2{t})(i) = mean(mean(PSTH_Inc_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t})/binSize,2)-mean(PSTH_Inc_All{1, VAR_SDixOIMatched(i)}(:,1:10)/binSize,2));
            CorFR_AUC.All.(trialTypo2{t})(i,:) = sum(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t})/0.2,1)-mean(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,1:10)/0.2,1)));
            IncFR_AUC.All.(trialTypo2{t})(i,:) = sum(mean(PSTH_Inc_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t})/0.2,1)-mean(mean(PSTH_Inc_All{1, VAR_SDixOIMatched(i)}(:,1:10)/0.2,1)));
    end
end

%% 7. Get cue onset type distribution
for ordType = [2 3 4]
    if ordType == 2 | ordType == 3
        l = ordType-1;
    elseif ordType == 4
        l = 3;
    end
    SD_All_clust.(trialTypo2{l}) = All{28, 1}(neuronIndex.PyrIxVAR_SD,ordType);
    SD_OI_clust.(trialTypo2{l}) = All{28, 1}(neuronIndex.PyrIx_VAR_SD_OI,ordType);
end

[B_OI, ~, ib_OI] = unique(All{28, 1}(neuronIndex.PyrIx_VAR_SD_OI,[2 3 4]), 'rows');
numoccurences_OI = accumarray(ib_OI, 1);

dIx = [];
upIx = [];
for i = 1:size(B_OI,1)
    if numel(find(B_OI(i,2:3)==[-1 -1]))==2
        dIx = [dIx i];
    elseif numel(find(B_OI(i,2:3)==[1 1]))==2
        upIx = [upIx i];
    end
end
for i = 1:size(neuronIndex.PyrIx_VAR_SD_OI,2)
    if ismember(ib_OI(i),dIx)
        ClustType_OI(i) = -1;
    elseif ismember(ib_OI(i),upIx)
        ClustType_OI(i) = 1;
    else
        ClustType_OI(i) = 0;
    end
end

[B_All, ~, ib_All] = unique(All{28, 1}(neuronIndex.PyrIxVAR_SD,[2 3 4]), 'rows');
numoccurences_All = accumarray(ib_All, 1);
dIx = [];
upIx = [];
for i = 1:size(B_All,1)
    if numel(find(B_All(i,2:3)==[-1 -1]))==2
        dIx = [dIx i];
    elseif numel(find(B_All(i,2:3)==[1 1]))==2
        upIx = [upIx i];
    end
end
for i = 1:size(neuronIndex.PyrIxVAR_SD,2)
    if ismember(ib_All(i),dIx)
        ClustType_All(i) = -1;
    elseif ismember(ib_All(i),upIx)
        ClustType_All(i) = 1;
    else
        ClustType_All(i) = 0;
    end
end

SD_OI_clust.TS = ClustType_OI';
SD_OI_clust.TH = ClustType_OI';
SD_OI_clust.cue_bef = ClustType_OI';

SD_All_clust.TS = ClustType_All';
SD_All_clust.TH = ClustType_All';
SD_All_clust.cue_bef = ClustType_All';

ClustUpMidDown = [1 0 -1];
tot_All_OI = numel(SD_OI_clust.TS);
tot_All_All = numel(VAR_SDixOIMatched);
for i = 1:numel(ClustUpMidDown)
    x = [zeros(numel(SD_OI_clust.TS),1); ones(numel(VAR_SDixOIMatched),1)]';
    y = [zeros(numel(find(SD_OI_clust.TS == ClustUpMidDown(i))),1); ones(tot_All_OI-numel(find(SD_OI_clust.TS == ClustUpMidDown(i))),1); zeros(numel(find(SD_All_clust.TS(VAR_SDixOIMatched) == ClustUpMidDown(i))),1); ones(tot_All_All-numel(find(SD_All_clust.TS(VAR_SDixOIMatched) == ClustUpMidDown(i))),1)]';
    [tbl_cueTypesOIvsOther{i},chi2stat_cueTypesOIvsOther{i},Chi_cueTypesOIvsOther.pval(i)] = crosstab(x,y);
    if numel(find(SD_OI_clust.TS == ClustUpMidDown(i))) < 5 || numel(find(SD_All_clust.TS(VAR_SDixOIMatched) == ClustUpMidDown(i))) < 5
        tbl_cueTypesOIvsOther{i} = table([numel(find(SD_OI_clust.TS == ClustUpMidDown(i)));tot_All_OI-numel(find(SD_OI_clust.TS == ClustUpMidDown(i)))],[numel(find(SD_All_clust.TS(VAR_SDixOIMatched) == ClustUpMidDown(i)));tot_All_All-numel(find(SD_All_clust.TS(VAR_SDixOIMatched) == ClustUpMidDown(i)))],'VariableNames',{'Flu','NoFlu'},'RowNames',{'NoShot','Shot'});
        [~,Chi_cueTypesOIvsOther.pval(i),chi2stat_cueTypesOIvsOther{i}] = fishertest(tbl_cueTypesOIvsOther{i});
    end
end

%% 3. Get boxplots with signrank comp per sync point and up, down or unm
n = 1;
pos = [1 2 3];
for t = [1 2 3]
        for Cl = 1:numel(ClustUpMidDown)
            diffCorInc.OI.(trialTypo2{t})(Cl) = median(((IncFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))'-CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))')./abs((CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))')))*100);

            try
            [CorVsIncSign.OI.(trialTypo2{t})(Cl), ~, CorVsIncSign.OIstats.(trialTypo2{t})(Cl)] = signrank(CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),IncFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))));
            catch
                CorVsIncSign.OI.(trialTypo2{t})(Cl) = nan;
            end
        end
        n = n+1;
        [CorVsIncSign_Combined.OI.(trialTypo2{t}), ~, CorVsIncSign_Combined.OI.stats.(trialTypo2{t})] = signrank(abs(CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))),abs(IncFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))));
        CorVsIncSign_Combined2.OI.(trialTypo2{t}){2} = abs(CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2))));
        CorVsIncSign_Combined2.OI.(trialTypo2{t}){3} = abs(IncFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2))));
end
CorVsIncSign_Combined.OI_combined = mafdr([CorVsIncSign_Combined.OI.TH CorVsIncSign_Combined.OI.Cue],'BHFDR',true);
CorVsIncSign_Combined.OI.TH = CorVsIncSign_Combined.OI_combined(1);
CorVsIncSign_Combined.OI.Cue = CorVsIncSign_Combined.OI_combined(2);

figure
subplot(1,3,1)
hold on
plot([abs(CorFR.OI.TS(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1)))); abs(IncFR.OI.TS(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1))))],'r');
plot([abs(CorFR.OI.TS(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3)))); abs(IncFR.OI.TS(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3))))],'b');
boxplot([CorVsIncSign_Combined2.OI.TS{1, 2}' CorVsIncSign_Combined2.OI.TS{1, 3}'],'Labels',{'Correct','Incature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 20])
subplot(1,3,2)
hold on
plot([abs(CorFR.OI.TH(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1)))); abs(IncFR.OI.TH(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1))))],'r');
plot([abs(CorFR.OI.TH(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3)))); abs(IncFR.OI.TH(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3))))],'b');
boxplot([CorVsIncSign_Combined2.OI.TH{1, 2}' CorVsIncSign_Combined2.OI.TH{1, 3}'],'Labels',{'Correct','Incature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 20])
subplot(1,3,3)
hold on
plot([abs(CorFR.OI.Cue(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1)))); abs(IncFR.OI.Cue(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1))))],'r');
plot([abs(CorFR.OI.Cue(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3)))); abs(IncFR.OI.Cue(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3))))],'b');
boxplot([CorVsIncSign_Combined2.OI.Cue{1, 2}' CorVsIncSign_Combined2.OI.Cue{1, 3}'],'Labels',{'Correct','Incature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 20])
