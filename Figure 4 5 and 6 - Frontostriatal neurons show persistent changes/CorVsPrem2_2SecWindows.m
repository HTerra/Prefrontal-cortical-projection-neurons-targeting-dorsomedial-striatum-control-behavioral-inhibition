function CorVsPrem2_2SecWindows(All,neuronIndex, SD_All, ITI_All, SD_OI, ITI_OI, ITIperiod)

% 1st figure = Figure 6C
% 4th figure = Figure 7C
% 7th figure = Figure 4H
% 8th figure = Figure 4G
% 9th figure = Figure 4I
% 10th figure = Figure 6B
% 11th figure = Figure 6A

binSize = 0.2;
ITIperiod = 'ITI125';
trialTypo = {'TS','TH','Resp'};
trialTypo2 = {'TS','TH','Resp_bef'};
timeBins = {11:15, 38:47, 51:60};

%% 1. Get CorFr for OI and All per sync point

for i = 1:numel(neuronIndex.PyrIx_VAR_ITI_OI)
    PSTH_Cor_OI{i} = [ITI_OI.(ITIperiod).psthBinsValue{1, i}.Cor.trial_start(:,1:25) ITI_OI.(ITIperiod).psthBinsValue{1, i}.Cor.wait_start(:,1:25) ITI_OI.(ITIperiod).psthBinsValue{1, i}.Cor.resp(:,1:20)];
    PSTH_Prem_OI{i} = [ITI_OI.(ITIperiod).psthBinsValue{1, i}.Prem.trial_start(:,1:25) ITI_OI.(ITIperiod).psthBinsValue{1, i}.Prem.wait_start(:,1:25) ITI_OI.(ITIperiod).psthBinsValue{1, i}.Prem.resp(:,1:20)];
    TresholdTime_Cor_OI(i) = mean(ITI_OI.(ITIperiod).TrialDist{1, i}.CorTreshold);
    TresholdTime_Prem_OI(i) = mean(ITI_OI.(ITIperiod).TrialDist{1, i}.PremTreshold);
    TresholdTime.Cor_Raw_OI{i} = ITI_OI.(ITIperiod).TrialDist{1, i}.CorTreshold;
    TresholdTime.Prem_Raw_OI{i} = ITI_OI.(ITIperiod).TrialDist{1, i}.PremTreshold;
end

for i = 1:numel(neuronIndex.PyrIxVAR_ITI)
    PSTH_Cor_All{i} = [ITI_All.(ITIperiod).psthBinsValue{1, i}.Cor.trial_start(:,1:25) ITI_All.(ITIperiod).psthBinsValue{1, i}.Cor.wait_start(:,1:25) ITI_All.(ITIperiod).psthBinsValue{1, i}.Cor.resp(:,1:20)];
    PSTH_Prem_All{i} = [ITI_All.(ITIperiod).psthBinsValue{1, i}.Prem.trial_start(:,1:25) ITI_All.(ITIperiod).psthBinsValue{1, i}.Prem.wait_start(:,1:25) ITI_All.(ITIperiod).psthBinsValue{1, i}.Prem.resp(:,1:20)];
    TresholdTime_Cor_All(i) = mean(ITI_All.(ITIperiod).TrialDist{1, i}.CorTreshold);
    TresholdTime_Prem_All(i) = mean(ITI_All.(ITIperiod).TrialDist{1, i}.PremTreshold);
    TresholdTime.Cor_Raw_All{i} = ITI_All.(ITIperiod).TrialDist{1, i}.CorTreshold;
    TresholdTime.Prem_Raw_All{i} = ITI_All.(ITIperiod).TrialDist{1, i}.PremTreshold;
end

%% 1. Get CorFr for OI and All per sync point, calculate significance

for t = 1:size(timeBins,2)
    for i = 1:numel(neuronIndex.PyrIx_VAR_ITI_OI)

                p = ranksum(mean(PSTH_Cor_OI{1, i}(:,timeBins{t}),2)/binSize,mean(PSTH_Prem_OI{1, i}(:,timeBins{t})/binSize,2));
                Sign.ITI_OI.(trialTypo2{t})(i) = p;

            CorFR.OI.(trialTypo2{t})(i) = mean(mean(PSTH_Cor_OI{1, i}(:,timeBins{t})/binSize,2)-mean(PSTH_Cor_OI{1, i}(:,1:10)/binSize,2));
            PremFR.OI.(trialTypo2{t})(i) = mean(mean(PSTH_Prem_OI{1, i}(:,timeBins{t})/binSize,2)-mean(PSTH_Prem_OI{1, i}(:,1:10)/binSize,2));
            CorFR_AUC.OI.(trialTypo2{t})(i,:) = sum(mean(PSTH_Cor_OI{1, i}(:,timeBins{t})/binSize,1)-mean(mean(PSTH_Cor_OI{1, i}(:,1:10)/binSize,1)));
            PremFR_AUC.OI.(trialTypo2{t})(i,:) = sum(mean(PSTH_Prem_OI{1, i}(:,timeBins{t})/binSize,1)-mean(mean(PSTH_Prem_OI{1, i}(:,1:10)/binSize,1)));
    end
end

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

for t = 1:size(timeBins,2)
    for i = 1:numel(neuronIndex.PyrIxVAR_ITI(VAR_ITIixOIMatched))
                p = ranksum(mean(PSTH_Cor_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t}),2)/binSize,mean(PSTH_Prem_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t}),2)/binSize);
                Sign.ITI_All.(trialTypo2{t})(i) = p;
            CorFR.All.(trialTypo2{t})(i) = mean(mean(PSTH_Cor_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t})/binSize,2)-mean(PSTH_Cor_All{1, VAR_ITIixOIMatched(i)}(:,1:10)/binSize,2));
            PremFR.All.(trialTypo2{t})(i) = mean(mean(PSTH_Prem_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t})/binSize,2)-mean(PSTH_Prem_All{1, VAR_ITIixOIMatched(i)}(:,1:10)/binSize,2));
            CorFR_AUC.All.(trialTypo2{t})(i,:) = sum(mean(PSTH_Cor_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t})/binSize,1)-mean(mean(PSTH_Cor_All{1, VAR_ITIixOIMatched(i)}(:,1:10)/binSize,1)));
            PremFR_AUC.All.(trialTypo2{t})(i,:) = sum(mean(PSTH_Prem_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t})/binSize,1)-mean(mean(PSTH_Prem_All{1, VAR_ITIixOIMatched(i)}(:,1:10)/binSize,1)));
    end
end

%% 7. Get Response type distribution
for ordType = [2 3 4]
    if ordType == 2 | ordType == 3
        l = ordType-1;
    elseif ordType == 4
        l = 3;
    end
    ITI_All_clust.(trialTypo2{l}) = All{28, 1}(neuronIndex.PyrIxVAR_ITI,ordType);
    ITI_OI_clust.(trialTypo2{l}) = All{28, 1}(neuronIndex.PyrIx_VAR_ITI_OI,ordType);
end

[B_OI, ~, ib_OI] = unique(All{28, 1}(neuronIndex.PyrIx_VAR_ITI_OI,[2 3 4]), 'rows');
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
for i = 1:size(neuronIndex.PyrIx_VAR_ITI_OI,2)
    if ismember(ib_OI(i),dIx)
        ClustType_OI(i) = -1;
    elseif ismember(ib_OI(i),upIx)
        ClustType_OI(i) = 1;
    else
        ClustType_OI(i) = 0;
    end
end

[B_All, ~, ib_All] = unique(All{28, 1}(neuronIndex.PyrIxVAR_ITI,[2 3 4]), 'rows');
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
for i = 1:size(neuronIndex.PyrIxVAR_ITI,2)
    if ismember(ib_All(i),dIx)
        ClustType_All(i) = -1;
    elseif ismember(ib_All(i),upIx)
        ClustType_All(i) = 1;
    else
        ClustType_All(i) = 0;
    end
end

ITI_OI_clust.TS = ClustType_OI';
ITI_OI_clust.TH = ClustType_OI';
ITI_OI_clust.Resp_bef = ClustType_OI';

ITI_All_clust.TS = ClustType_All';
ITI_All_clust.TH = ClustType_All';
ITI_All_clust.Resp_bef = ClustType_All';

%% Figure 6C
figure
bar([numel(find(ITI_All_clust.TS(VAR_ITIixOIMatched,1)==1))/numel(VAR_ITIixOIMatched) numel(find(ITI_OI_clust.TS==1))/numel(ITI_OI_clust.TS); numel(find(ITI_All_clust.TS(VAR_ITIixOIMatched,1)==0))/numel(VAR_ITIixOIMatched) numel(find(ITI_OI_clust.TS==0))/numel(ITI_OI_clust.TS); numel(find(ITI_All_clust.TS(VAR_ITIixOIMatched,1)==-1))/numel(VAR_ITIixOIMatched) numel(find(ITI_OI_clust.TS==-1))/numel(ITI_OI_clust.TS)]); 
legend({'Other', 'Frontostriatal'})

ClustUpMidDown = [1 0 -1];
tot_All_OI = numel(ITI_OI_clust.TS);
tot_All_All = numel(VAR_ITIixOIMatched);
for i = 1:numel(ClustUpMidDown)
    x = [zeros(numel(ITI_OI_clust.TS),1); ones(numel(VAR_ITIixOIMatched),1)]';
    y = [zeros(numel(find(ITI_OI_clust.TS == ClustUpMidDown(i))),1); ones(tot_All_OI-numel(find(ITI_OI_clust.TS == ClustUpMidDown(i))),1); zeros(numel(find(ITI_All_clust.TS(VAR_ITIixOIMatched) == ClustUpMidDown(i))),1); ones(tot_All_All-numel(find(ITI_All_clust.TS(VAR_ITIixOIMatched) == ClustUpMidDown(i))),1)]';
    [tbl_RespTypesOIvsOther{i},chi2stat_RespTypesOIvsOther{i},Chi_RespTypesOIvsOther.pval(i)] = crosstab(x,y);
    if numel(find(ITI_OI_clust.TS == ClustUpMidDown(i))) < 5 || numel(find(ITI_All_clust.TS(VAR_ITIixOIMatched) == ClustUpMidDown(i))) < 5
        tbl_RespTypesOIvsOther{i} = table([numel(find(ITI_OI_clust.TS == ClustUpMidDown(i)));tot_All_OI-numel(find(ITI_OI_clust.TS == ClustUpMidDown(i)))],[numel(find(ITI_All_clust.TS(VAR_ITIixOIMatched) == ClustUpMidDown(i)));tot_All_All-numel(find(ITI_All_clust.TS(VAR_ITIixOIMatched) == ClustUpMidDown(i)))],'VariableNames',{'Flu','NoFlu'},'RowNames',{'NoShot','Shot'});
        [~,Chi_RespTypesOIvsOther.pval(i),chi2stat_RespTypesOIvsOther{i}] = fishertest(tbl_RespTypesOIvsOther{i});
    end
end

%% 2. Make scatterplots per sync point per up, down and unm

ClustUpMidDown = [1 0 -1];
figure
n = 1;
pos = [1 4 7 2 5 8 3 6 9];
for t = [1 2 3]
    for Cl = 1:numel(ClustUpMidDown)
        subplot(3,3,pos(n))
        hold on
        scatter(CorFR.OI.(trialTypo2{t})(find(Sign.ITI_OI.(trialTypo2{t})(1,:)>=0.05&ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),PremFR.OI.(trialTypo2{t})(find(Sign.ITI_OI.(trialTypo2{t})(1,:)>=0.05&ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),'r')
        scatter(CorFR.OI.(trialTypo2{t})(find(Sign.ITI_OI.(trialTypo2{t})(1,:)<0.05&ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),PremFR.OI.(trialTypo2{t})(find(Sign.ITI_OI.(trialTypo2{t})(1,:)<0.05&ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),'r','filled')
        %lsline
        scatter(CorFR.All.(trialTypo2{t})(find(Sign.ITI_All.(trialTypo2{t})(1,:)<0.05&ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))),PremFR.All.(trialTypo2{t})(find(Sign.ITI_All.(trialTypo2{t})(1,:)<0.05&ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))),'b','filled')
        scatter(CorFR.All.(trialTypo2{t})(find(Sign.ITI_All.(trialTypo2{t})(1,:)>=0.05&ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))),PremFR.All.(trialTypo2{t})(find(Sign.ITI_All.(trialTypo2{t})(1,:)>=0.05&ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))),'b')
        m = max([CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))) PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))) CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))) PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))]);
        m2 = min([CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))) PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))) CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))) PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))]);
        [R P] = corrcoef(CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))));
        R_OI(n) = R(2);
        R_P_OI(n) = P(2);
        [R P] = corrcoef(CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))),PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))));
        R_All(n) = R(2);
        R_P_All(n) = P(2);
        line([m2 m],[m2 m])
        xlim([m2 m])
        ylim([m2 m])
        xlabel('Cor (Hz)')
        ylabel('Prem (Hz)')
        n = n+1;
    end
end

%% 3. Get boxplots with signrank comp per sync point and up, down or unm
figure
n = 1;
pos = [1 2 3];
for t = [1 2 3]
        for Cl = 1:numel(ClustUpMidDown)
            hold on
            diffCorPrem.OI.(trialTypo2{t})(Cl) = median(((PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))'-CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))')./abs((CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))')))*100);

            try
            [CorVsPremSign.OI.(trialTypo2{t})(Cl), ~, CorVsPremSign.OIstats.(trialTypo2{t})(Cl)] = signrank(CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))));
            catch
                CorVsPremSign.OI.(trialTypo2{t})(Cl) = nan;
            end
        end
        subplot(3,1,pos(n))
        bar(diffCorPrem.OI.(trialTypo2{t}))
        ylim([-25 25])
        n = n+1;
        [CorVsPremSign_Combined.OI.(trialTypo2{t}), ~, CorVsPremSign_Combined.OI.stats.(trialTypo2{t})] = signrank(abs(CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))),abs(PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))));
        CorVsPremSign_Combined2.OI.(trialTypo2{t}){2} = abs(CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2))));
        CorVsPremSign_Combined2.OI.(trialTypo2{t}){3} = abs(PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2))));
end
CorVsPremSign_Combined.OI_combined = [CorVsPremSign_Combined.OI.TH CorVsPremSign_Combined.OI.Resp_bef];
CorVsPremSign_Combined.OI_combined = mafdr([CorVsPremSign_Combined.OI.TH CorVsPremSign_Combined.OI.Resp_bef],'BHFDR',true);
CorVsPremSign_Combined.OI.TH = CorVsPremSign_Combined.OI_combined(1);
CorVsPremSign_Combined.OI.Resp_bef = CorVsPremSign_Combined.OI_combined(2);

%% Figure 7C
figure
subplot(1,3,1)
hold on
plot([CorVsPremSign_Combined2.OI.TS{1, 2}(find(ITI_OI_clust.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))'==ClustUpMidDown(1))); CorVsPremSign_Combined2.OI.TS{1, 3}(find(ITI_OI_clust.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))'==ClustUpMidDown(1)))],'r')
plot([CorVsPremSign_Combined2.OI.TS{1, 2}(find(ITI_OI_clust.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))'==ClustUpMidDown(3))); CorVsPremSign_Combined2.OI.TS{1, 3}(find(ITI_OI_clust.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))'==ClustUpMidDown(3)))],'b')
boxplot([CorVsPremSign_Combined2.OI.TS{1, 2}' CorVsPremSign_Combined2.OI.TS{1, 3}'],'Labels',{'Correct','Premature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 25])
subplot(1,3,2)
hold on
plot([CorVsPremSign_Combined2.OI.TH{1, 2}(find(ITI_OI_clust.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))'==ClustUpMidDown(1))); CorVsPremSign_Combined2.OI.TH{1, 3}(find(ITI_OI_clust.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))'==ClustUpMidDown(1)))],'r')
plot([CorVsPremSign_Combined2.OI.TH{1, 2}(find(ITI_OI_clust.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))'==ClustUpMidDown(3))); CorVsPremSign_Combined2.OI.TH{1, 3}(find(ITI_OI_clust.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))'==ClustUpMidDown(3)))],'b')
boxplot([CorVsPremSign_Combined2.OI.TH{1, 2}' CorVsPremSign_Combined2.OI.TH{1, 3}'],'Labels',{'Correct','Premature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 25])
subplot(1,3,3)
hold on
plot([CorVsPremSign_Combined2.OI.Resp_bef{1, 2}(find(ITI_OI_clust.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))'==ClustUpMidDown(1))); CorVsPremSign_Combined2.OI.Resp_bef{1, 3}(find(ITI_OI_clust.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))'==ClustUpMidDown(1)))],'r')
plot([CorVsPremSign_Combined2.OI.Resp_bef{1, 2}(find(ITI_OI_clust.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))'==ClustUpMidDown(3))); CorVsPremSign_Combined2.OI.Resp_bef{1, 3}(find(ITI_OI_clust.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))'==ClustUpMidDown(3)))],'b')
boxplot([CorVsPremSign_Combined2.OI.Resp_bef{1, 2}' CorVsPremSign_Combined2.OI.Resp_bef{1, 3}'],'Labels',{'Correct','Premature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 25])

%%
figure
n = 1;
pos = [1 2 3];
for t = [1 2 3]
        for Cl = 1:numel(ClustUpMidDown)
            diffCorPrem.All.(trialTypo2{t})(Cl) = median(((PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))'-CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))')./abs((CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))')))*100);

            h1(i) = lillietest(CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))));
            h2(i) = lillietest(PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))); 
            if h1(i) == 0 & h2(i) == 0
                try
                [h CorVsPremSign.All.(trialTypo2{t})(Cl)] = ttest(CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))),PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))),'Vartype','unequal');
                catch
                    CorVsPremSign.All.(trialTypo2{t})(Cl) = nan;
                end
            else
                try
                CorVsPremSign.All.(trialTypo2{t})(Cl) = signrank(CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))),PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))));
                catch
                    CorVsPremSign.All.(trialTypo2{t})(Cl) = nan;
            end
        end
        end
        subplot(3,1,pos(n))
        bar(diffCorPrem.All.(trialTypo2{t}))
        ylim([-25 25])
        n = n+1;
        CorVsPremSign_Combined.All.(trialTypo2{t}) = signrank(abs(CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'~=ClustUpMidDown(2)))),abs(PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'~=ClustUpMidDown(2)))));
        CorVsPremSign_Combined2.All.(trialTypo2{t}){2} = abs(CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'~=ClustUpMidDown(2))));
        CorVsPremSign_Combined2.All.(trialTypo2{t}){3} = abs(PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'~=ClustUpMidDown(2))));
end
CorVsPremSign_Combined.All_combined = mafdr([CorVsPremSign_Combined.All.TH CorVsPremSign_Combined.All.Resp_bef],'BHFDR',true);

%%
figure
subplot(1,3,1)
hold on
plot([CorVsPremSign_Combined2.All.TS{1, 2}; CorVsPremSign_Combined2.All.TS{1, 3}],'k')
boxplot([CorVsPremSign_Combined2.All.TS{1, 2}' CorVsPremSign_Combined2.All.TS{1, 3}'],'Labels',{'Correct','Premature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 25])
subplot(1,3,2)
hold on
plot([CorVsPremSign_Combined2.All.TH{1, 2}; CorVsPremSign_Combined2.All.TH{1, 3}],'k')
boxplot([CorVsPremSign_Combined2.All.TH{1, 2}' CorVsPremSign_Combined2.All.TH{1, 3}'],'Labels',{'Correct','Premature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 25])
subplot(1,3,3)
hold on
plot([CorVsPremSign_Combined2.All.Resp_bef{1, 2}; CorVsPremSign_Combined2.All.Resp_bef{1, 3}],'k')
boxplot([CorVsPremSign_Combined2.All.Resp_bef{1, 2}' CorVsPremSign_Combined2.All.Resp_bef{1, 3}'],'Labels',{'Correct','Premature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 25])

% Get baseline subtracted PETHs
for n = 1:size(ITI_OI.(ITIperiod).psthBinsValue,2)
    OI_PSTH.PSTH.TS.Cor(n,:) = (mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:85)/binSize,1)-mean(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));% /max(abs(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,11:15),1)));
    OI_PSTH.PSTH.TH.Cor(n,:) = (mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.wait_start(:,1:25)/binSize,1)-mean(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));%/max(abs(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.wait_start(:,6:25),1)));
    OI_PSTH.PSTH.Resp_bef.Cor(n,:) = (mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.resp(:,1:30)/binSize,1)-mean(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));%/max(abs(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.resp(:,1:10),1)));
    OI_PSTH.PSTH.TS.Prem(n,:) = (mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Prem.trial_start(:,1:85)/binSize,1)-mean(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Prem.trial_start(:,1:10)/binSize,1));%/max(abs(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,11:15),1)));
    OI_PSTH.PSTH.TH.Prem(n,:) = (mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Prem.wait_start(:,1:25)/binSize,1)-mean(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Prem.trial_start(:,1:10)/binSize,1));%/max(abs(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.wait_start(:,6:25),1)));
    OI_PSTH.PSTH.Resp_bef.Prem(n,:) = (mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Prem.resp(:,1:30)/binSize,1)-mean(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Prem.trial_start(:,1:10)/binSize,1));%/max(abs(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Prem.resp(:,1:10),1)));
end

for n = 1:size(ITI_All.(ITIperiod).psthBinsValue,2)
    All_PSTH.PSTH.TS.Cor(n,:) = (mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:85)/binSize,1)-mean(mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));% /max(abs(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,11:15),1)));
    All_PSTH.PSTH.TH.Cor(n,:) = (mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Cor.wait_start(:,1:25)/binSize,1)-mean(mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));%/max(abs(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.wait_start(:,6:25),1)));
    All_PSTH.PSTH.Resp_bef.Cor(n,:) = (mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Cor.resp(:,1:30)/binSize,1)-mean(mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));%/max(abs(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.resp(:,1:10),1)));
    All_PSTH.PSTH.TS.Prem(n,:) = (mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Prem.trial_start(:,1:85)/binSize,1)-mean(mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Prem.trial_start(:,1:10)/binSize,1));%/max(abs(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,11:15),1)));
    All_PSTH.PSTH.TH.Prem(n,:) = (mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Prem.wait_start(:,1:25)/binSize,1)-mean(mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Prem.trial_start(:,1:10)/binSize,1));%/max(abs(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Cor.wait_start(:,6:25),1)));
    All_PSTH.PSTH.Resp_bef.Prem(n,:) = (mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Prem.resp(:,1:30)/binSize,1)-mean(mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(ITI_All.(ITIperiod).psthBinsValue{1, n}.Prem.trial_start(:,1:10)/binSize,1));%/max(abs(mean(ITI_OI.(ITIperiod).psthBinsValue{1, n}.Prem.resp(:,1:10),1)));
end

time_axis_TS = linspace(-2,15,((2+15)/binSize));
time_axis_TH = linspace(-2,3,((2+3)/binSize));
time_axis_Resp_bef = linspace(-2,2,((2+2)/binSize));
axisRange = [-6 6];

%[B,I] = sort(mean(OI_PSTH.PSTH.Resp_bef.Cor(:,1:5),2),'descend');
[B,I_one] = sort(mean(OI_PSTH.PSTH.Resp_bef.Cor(find(ClustType_OI==1),1:5),2),'descend');
[B,I_zero] = sort(mean(OI_PSTH.PSTH.Resp_bef.Cor(find(ClustType_OI==0),1:5),2),'descend');
[B,I_neg] = sort(mean(OI_PSTH.PSTH.Resp_bef.Cor(find(ClustType_OI==-1),1:5),2),'descend');
I_oneIxs = find(ClustType_OI==1);
I_zeroIxs = find(ClustType_OI==0);
I_negIxs = find(ClustType_OI==-1);
I = [I_oneIxs(I_one) I_zeroIxs(I_zero) I_negIxs(I_neg)];

%% Figure 4H
figure
colormap('jet')
imagesc(B_OI(ib_OI(I),1:3))

%% Figure 4G
figure
subplot(1,3,1)
colormap('jet')
%[B,I] = sort(mean(OI_PSTH.PSTH.TS.Cor(:,11:15),2),'descend');
imagesc(time_axis_TS,1:(size(OI_PSTH.PSTH.TS.Cor,1)),OI_PSTH.PSTH.TS.Cor(I,1:85))
axis([-2-(0.5*binSize) 15+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.TS.Cor,1)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Trial start')

subplot(1,3,2)
colormap('jet')
%[B,I] = sort(mean(OI_PSTH.PSTH.TH.Cor(:,11:15),2),'descend');
imagesc(time_axis_TH,1:(size(OI_PSTH.PSTH.TH.Cor,1)),OI_PSTH.PSTH.TH.Cor(I,1:25))
axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.TH.Cor,1)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Cue orientation')

subplot(1,3,3)
colormap('jet')
%[B,I] = sort(mean(OI_PSTH.PSTH.Resp_bef.Cor(:,1:10),2),'descend');
imagesc(time_axis_Resp_bef,1:(size(OI_PSTH.PSTH.Resp_bef.Cor,1)),OI_PSTH.PSTH.Resp_bef.Cor(I,1:20))
axis([-2-(0.5*binSize) 2+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.Resp_bef.Cor,1)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Response')

%% Figure 4I
figure
subplot(1,3,1) 
hold on
SEM = nanstd(OI_PSTH.PSTH.TS.Cor(I_oneIxs(I_one),1:85),1)/sqrt(size(OI_PSTH.PSTH.TS.Cor(I_oneIxs(I_one),1:85),1)); 
shadedErrorBar(time_axis_TS,mean(OI_PSTH.PSTH.TS.Cor(I_oneIxs(I_one),1:85)),SEM,{'r'});
SEM = nanstd(OI_PSTH.PSTH.TS.Cor(I_zeroIxs(I_zero),1:85),1)/sqrt(size(OI_PSTH.PSTH.TS.Cor(I_zeroIxs(I_zero),1:85),1)); 
shadedErrorBar(time_axis_TS,mean(OI_PSTH.PSTH.TS.Cor(I_zeroIxs(I_zero),1:85)),SEM,{'k'});
SEM = nanstd(OI_PSTH.PSTH.TS.Cor(I_negIxs(I_neg),1:85),1)/sqrt(size(OI_PSTH.PSTH.TS.Cor(I_negIxs(I_neg),1:85),1)); 
shadedErrorBar(time_axis_TS,mean(OI_PSTH.PSTH.TS.Cor(I_negIxs(I_neg),1:85)),SEM,{'b'});
subplot(1,3,2) 
hold on
SEM = nanstd(OI_PSTH.PSTH.TH.Cor(I_oneIxs(I_one),1:25),1)/sqrt(size(OI_PSTH.PSTH.TH.Cor(I_oneIxs(I_one),1:25),1)); 
shadedErrorBar(time_axis_TH,mean(OI_PSTH.PSTH.TH.Cor(I_oneIxs(I_one),1:25)),SEM,{'r'});
SEM = nanstd(OI_PSTH.PSTH.TH.Cor(I_zeroIxs(I_zero),1:25),1)/sqrt(size(OI_PSTH.PSTH.TH.Cor(I_zeroIxs(I_zero),1:25),1)); 
shadedErrorBar(time_axis_TH,mean(OI_PSTH.PSTH.TH.Cor(I_zeroIxs(I_zero),1:25)),SEM,{'k'});
SEM = nanstd(OI_PSTH.PSTH.TH.Cor(I_negIxs(I_neg),1:25),1)/sqrt(size(OI_PSTH.PSTH.TH.Cor(I_negIxs(I_neg),1:25),1)); 
shadedErrorBar(time_axis_TH,mean(OI_PSTH.PSTH.TH.Cor(I_negIxs(I_neg),1:25)),SEM,{'b'});
subplot(1,3,3) 
hold on
SEM = nanstd(OI_PSTH.PSTH.Resp_bef.Cor(I_oneIxs(I_one),1:20),1)/sqrt(size(OI_PSTH.PSTH.Resp_bef.Cor(I_oneIxs(I_one),1:20),1)); 
shadedErrorBar(time_axis_Resp_bef,mean(OI_PSTH.PSTH.Resp_bef.Cor(I_oneIxs(I_one),1:20)),SEM,{'r'});
SEM = nanstd(OI_PSTH.PSTH.Resp_bef.Cor(I_zeroIxs(I_zero),1:20),1)/sqrt(size(OI_PSTH.PSTH.Resp_bef.Cor(I_zeroIxs(I_zero),1:20),1)); 
shadedErrorBar(time_axis_Resp_bef,mean(OI_PSTH.PSTH.Resp_bef.Cor(I_zeroIxs(I_zero),1:20)),SEM,{'k'});
SEM = nanstd(OI_PSTH.PSTH.Resp_bef.Cor(I_negIxs(I_neg),1:20),1)/sqrt(size(OI_PSTH.PSTH.Resp_bef.Cor(I_negIxs(I_neg),1:20),1)); 
shadedErrorBar(time_axis_Resp_bef,mean(OI_PSTH.PSTH.Resp_bef.Cor(I_negIxs(I_neg),1:20)),SEM,{'b'});

I_oneIxs = find(ClustType_All(VAR_ITIixOIMatched)==1);
I_zeroIxs = find(ClustType_All(VAR_ITIixOIMatched)==0);
I_negIxs = find(ClustType_All(VAR_ITIixOIMatched)==-1);
[B,I_one] = sort(mean(All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched(I_oneIxs),1:5),2),'descend');
[B,I_zero] = sort(mean(All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched(I_zeroIxs),1:5),2),'descend');
[B,I_neg] = sort(mean(All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched(I_negIxs),1:5),2),'descend');
I = [I_oneIxs(I_one) I_zeroIxs(I_zero) I_negIxs(I_neg)];

%% Figure 6B
figure
colormap('jet')
imagesc(B_All(ib_All(VAR_ITIixOIMatched(I)),1:3))

%% Figure 6A
figure
subplot(1,3,1)
colormap('jet')
%[B,I] = sort(mean(All_PSTH.PSTH.TS.Cor(VAR_ITIixOIMatched,11:15),2),'descend');
imagesc(time_axis_TS,1:(size(All_PSTH.PSTH.TS.Cor(VAR_ITIixOIMatched),2)),All_PSTH.PSTH.TS.Cor(VAR_ITIixOIMatched(I),1:85))
axis([-2-(0.5*binSize) 15+(0.5*binSize) 0.5 size(All_PSTH.PSTH.TS.Cor(VAR_ITIixOIMatched),2)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Trial start')

subplot(1,3,2)
colormap('jet')
%[B,I] = sort(mean(All_PSTH.PSTH.TH.Cor(VAR_ITIixOIMatched,11:15),2),'descend');
imagesc(time_axis_TH,1:(size(All_PSTH.PSTH.TH.Cor(VAR_ITIixOIMatched),2)),All_PSTH.PSTH.TH.Cor(VAR_ITIixOIMatched(I),1:25))
axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(All_PSTH.PSTH.TH.Cor(VAR_ITIixOIMatched),2)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Cue orientation')

subplot(1,3,3)
colormap('jet')
%[B,I] = sort(mean(All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched,1:10),2),'descend');
imagesc(time_axis_Resp_bef,1:(size(All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched),2)),All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched(I),1:20))
axis([-2-(0.5*binSize) 2+(0.5*binSize) 0.5 size(All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched),2)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Response')
