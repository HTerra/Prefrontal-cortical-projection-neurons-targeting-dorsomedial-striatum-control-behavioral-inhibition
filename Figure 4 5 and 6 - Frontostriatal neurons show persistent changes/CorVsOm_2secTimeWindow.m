function CorVsOm_2secTimeWindow(All,neuronIndex, SD_All, SD_OI, SDperiod)

% Correct vs omission

% 1st figure = Figure S2D
% 4th figure = Figure S4C
% 7th figure = Figure 5E
% 8th figure = Figure 5D
% 9th figure = Figure 5F
% 10th figure = Figure S2B left
% 11th figure = Figure S2A

binSize = 0.2;
SDperiod = 'All';
trialTypo = {'TS','TH','Cue'};
trialTypo2 = {'TS','TH','Cue'};
timeBins = {11:15, 48:57, 26:35};

%% 1. Get CorFr for OI and All per sync point

for i = 1:numel(neuronIndex.PyrIx_VAR_SD_OI)
    PSTH_Cor_OI{i} = [SD_OI.(SDperiod).psthBinsValue{1, i}.Cor.trial_start(:,1:35) SD_OI.(SDperiod).psthBinsValue{1, i}.Cor.wait_start(:,1:25)];
    PSTH_Om_OI{i} = [SD_OI.(SDperiod).psthBinsValue{1, i}.Om.trial_start(:,1:35) SD_OI.(SDperiod).psthBinsValue{1, i}.Om.wait_start(:,1:25)];
end

for i = 1:numel(neuronIndex.PyrIxVAR_SD)
    PSTH_Cor_All{i} = [SD_All.(SDperiod).psthBinsValue{1, i}.Cor.trial_start(:,1:35) SD_All.(SDperiod).psthBinsValue{1, i}.Cor.wait_start(:,1:25)];
    PSTH_Om_All{i} = [SD_All.(SDperiod).psthBinsValue{1, i}.Om.trial_start(:,1:35) SD_All.(SDperiod).psthBinsValue{1, i}.Om.wait_start(:,1:25)];
end

%% 1. Get CorFr for OI and All per sync point, calculate significance

for t = 1:size(timeBins,2)
    for i = 1:numel(neuronIndex.PyrIx_VAR_SD_OI)
            p = ranksum(mean(PSTH_Cor_OI{1, i}(:,timeBins{t}),2)/binSize,mean(PSTH_Om_OI{1, i}(:,timeBins{t})/binSize,2));
            Sign.SD_OI.(trialTypo2{t})(i) = p;
            CorFR.OI.(trialTypo2{t})(i) = mean(mean(PSTH_Cor_OI{1, i}(:,timeBins{t})/binSize,2)-mean(PSTH_Cor_OI{1, i}(:,1:10)/binSize,2));
            OmFR.OI.(trialTypo2{t})(i) = mean(mean(PSTH_Om_OI{1, i}(:,timeBins{t})/binSize,2)-mean(PSTH_Om_OI{1, i}(:,1:10)/binSize,2));
            CorFR_AUC.OI.(trialTypo2{t})(i,:) = sum(mean(PSTH_Cor_OI{1, i}(:,timeBins{t})/binSize,1)-mean(mean(PSTH_Cor_OI{1, i}(:,1:10)/binSize,1)));
            OmFR_AUC.OI.(trialTypo2{t})(i,:) = sum(mean(PSTH_Om_OI{1, i}(:,timeBins{t})/binSize,1)-mean(mean(PSTH_Om_OI{1, i}(:,1:10)/binSize,1)));
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
            p = ranksum(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t}),2)/binSize,mean(PSTH_Om_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t}),2)/binSize);
            Sign.SD_All.(trialTypo2{t})(i) = p;
            CorFR.All.(trialTypo2{t})(i) = mean(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t})/binSize,2)-mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,1:10)/binSize,2));
            OmFR.All.(trialTypo2{t})(i) = mean(mean(PSTH_Om_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t})/binSize,2)-mean(PSTH_Om_All{1, VAR_SDixOIMatched(i)}(:,1:10)/binSize,2));
            CorFR_AUC.All.(trialTypo2{t})(i,:) = sum(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t})/0.2,1)-mean(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,1:10)/0.2,1)));
            OmFR_AUC.All.(trialTypo2{t})(i,:) = sum(mean(PSTH_Om_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t})/0.2,1)-mean(mean(PSTH_Om_All{1, VAR_SDixOIMatched(i)}(:,1:10)/0.2,1)));
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

%% Figure S2D

figure
bar([numel(find(SD_All_clust.TS(VAR_SDixOIMatched,1)==1))/numel(VAR_SDixOIMatched) numel(find(SD_OI_clust.TS==1))/numel(SD_OI_clust.TS); numel(find(SD_All_clust.TS(VAR_SDixOIMatched,1)==0))/numel(VAR_SDixOIMatched) numel(find(SD_OI_clust.TS==0))/numel(SD_OI_clust.TS); numel(find(SD_All_clust.TS(VAR_SDixOIMatched,1)==-1))/numel(VAR_SDixOIMatched) numel(find(SD_OI_clust.TS==-1))/numel(SD_OI_clust.TS)]); 
legend({'Other','Frontostriatal'})

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

%% 2. Make scatterplots per sync point per up, down and unm

ClustUpMidDown = [1 0 -1];
figure
n = 1;
pos = [1 4 7 2 5 8 3 6 9];
for t = [1 2 3]
    for Cl = 1:numel(ClustUpMidDown)
        subplot(3,3,pos(n))
        hold on
        scatter(CorFR.OI.(trialTypo2{t})(find(Sign.SD_OI.(trialTypo2{t})(1,:)>=0.05&SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),OmFR.OI.(trialTypo2{t})(find(Sign.SD_OI.(trialTypo2{t})(1,:)>=0.05&SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),'r')
        scatter(CorFR.OI.(trialTypo2{t})(find(Sign.SD_OI.(trialTypo2{t})(1,:)<0.05&SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),OmFR.OI.(trialTypo2{t})(find(Sign.SD_OI.(trialTypo2{t})(1,:)<0.05&SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),'r','filled')
        %lsline
        scatter(CorFR.All.(trialTypo2{t})(find(Sign.SD_All.(trialTypo2{t})(1,:)<0.05&SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))),OmFR.All.(trialTypo2{t})(find(Sign.SD_All.(trialTypo2{t})(1,:)<0.05&SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))),'b','filled')
        scatter(CorFR.All.(trialTypo2{t})(find(Sign.SD_All.(trialTypo2{t})(1,:)>=0.05&SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))),OmFR.All.(trialTypo2{t})(find(Sign.SD_All.(trialTypo2{t})(1,:)>=0.05&SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))),'b')
        m = max([CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))) OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))) CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))) OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))]);
        m2 = min([CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))) OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))) CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))) OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))]);
        [R P] = corrcoef(CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))));
        R_OI(n) = R(2);
        R_P_OI(n) = P(2);
        [R P] = corrcoef(CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))),OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))));
        R_All(n) = R(2);
        R_P_All(n) = P(2);
        line([m2 m],[m2 m])
        xlim([m2 m])
        ylim([m2 m])
        xlabel('Cor (Hz)')
        ylabel('Om (Hz)')
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
            diffCorOm.OI.(trialTypo2{t})(Cl) = median(((OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))'-CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))')./abs((CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))')))*100);
            try
            [CorVsOmSign.OI.(trialTypo2{t})(Cl), ~, CorVsOmSign.OIstats.(trialTypo2{t})(Cl)] = signrank(CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))));
            catch
                CorVsOmSign.OI.(trialTypo2{t})(Cl) = nan;
            end
        end
        subplot(3,1,pos(n))
        bar(diffCorOm.OI.(trialTypo2{t}))
        ylim([-25 25])
        n = n+1;
        [CorVsOmSign_Combined.OI.(trialTypo2{t}), ~, CorVsOmSign_Combined.OI.stats.(trialTypo2{t})] = signrank(abs(CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))),abs(OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))));
        CorVsOmSign_Combined2.OI.(trialTypo2{t}){2} = abs(CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2))));
        CorVsOmSign_Combined2.OI.(trialTypo2{t}){3} = abs(OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2))));
end
CorVsOmSign_Combined.OI_combined = mafdr([CorVsOmSign_Combined.OI.TH CorVsOmSign_Combined.OI.Cue],'BHFDR',true);;
CorVsOmSign_Combined.OI.TH = CorVsOmSign_Combined.OI_combined(1);
CorVsOmSign_Combined.OI.Cue = CorVsOmSign_Combined.OI_combined(2);

%% Figure S4D

figure
subplot(1,3,1)
hold on
plot([abs(CorFR.OI.TS(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1)))); abs(OmFR.OI.TS(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1))))],'r');
plot([abs(CorFR.OI.TS(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3)))); abs(OmFR.OI.TS(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3))))],'b');
boxplot([CorVsOmSign_Combined2.OI.TS{1, 2}' CorVsOmSign_Combined2.OI.TS{1, 3}'],'Labels',{'Correct','Omission'})
ylabel('Absolute delta FR (Hz)')
ylim([0 20])
subplot(1,3,2)
hold on
plot([abs(CorFR.OI.TH(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1)))); abs(OmFR.OI.TH(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1))))],'r');
plot([abs(CorFR.OI.TH(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3)))); abs(OmFR.OI.TH(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3))))],'b');
boxplot([CorVsOmSign_Combined2.OI.TH{1, 2}' CorVsOmSign_Combined2.OI.TH{1, 3}'],'Labels',{'Correct','Omission'})
ylabel('Absolute delta FR (Hz)')
ylim([0 20])
subplot(1,3,3)
hold on
plot([abs(CorFR.OI.Cue(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1)))); abs(OmFR.OI.Cue(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1))))],'r');
plot([abs(CorFR.OI.Cue(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3)))); abs(OmFR.OI.Cue(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3))))],'b');
boxplot([CorVsOmSign_Combined2.OI.Cue{1, 2}' CorVsOmSign_Combined2.OI.Cue{1, 3}'],'Labels',{'Correct','Omission'})
ylabel('Absolute delta FR (Hz)')
ylim([0 20])

%%

figure
n = 1;
pos = [1 2 3];
for t = [1 2 3]
        for Cl = 1:numel(ClustUpMidDown)
            diffCorOm.All.(trialTypo2{t})(Cl) = median(((OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))'-CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))')./abs((CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))')))*100);
            h1(i) = lillietest(CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))));
            h2(i) = lillietest(OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))); 
            if h1(i) == 0 & h2(i) == 0
                try
                [h CorVsOmSign.All.(trialTypo2{t})(Cl)] = ttest(CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))),OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))),'Vartype','unequal');
                catch
                    CorVsOmSign.All.(trialTypo2{t})(Cl) = nan;
                end
            else
                try
                CorVsOmSign.All.(trialTypo2{t})(Cl) = signrank(CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))),OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))));
                catch
                    CorVsOmSign.All.(trialTypo2{t})(Cl) = nan;
            end
        end
        end
        subplot(3,1,pos(n))
        bar(diffCorOm.All.(trialTypo2{t}))
        ylim([-25 25])
        n = n+1;
        CorVsOmSign.All.(trialTypo2{t}) = mafdr([CorVsOmSign.All.(trialTypo2{t})],'BHFDR',true);
        CorVsOmSign_Combined.All.(trialTypo2{t}) = signrank(abs(CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'~=ClustUpMidDown(2)))),abs(OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'~=ClustUpMidDown(2)))));
        CorVsOmSign_Combined2.All.(trialTypo2{t}){2} = abs(CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'~=ClustUpMidDown(2))));
        CorVsOmSign_Combined2.All.(trialTypo2{t}){3} = abs(OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'~=ClustUpMidDown(2))));
end

%%

figure
subplot(1,3,1)
hold on
plot([CorVsOmSign_Combined2.All.TS{1, 2}; CorVsOmSign_Combined2.All.TS{1, 3}],'k')
boxplot([CorVsOmSign_Combined2.All.TS{1, 2}' CorVsOmSign_Combined2.All.TS{1, 3}'],'Labels',{'Correct','Omature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 25])
subplot(1,3,2)
hold on
plot([CorVsOmSign_Combined2.All.TH{1, 2}; CorVsOmSign_Combined2.All.TH{1, 3}],'k')
boxplot([CorVsOmSign_Combined2.All.TH{1, 2}' CorVsOmSign_Combined2.All.TH{1, 3}'],'Labels',{'Correct','Omature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 25])
subplot(1,3,3)
hold on
plot([CorVsOmSign_Combined2.All.Cue{1, 2}; CorVsOmSign_Combined2.All.Cue{1, 3}],'k')
boxplot([CorVsOmSign_Combined2.All.Cue{1, 2}' CorVsOmSign_Combined2.All.Cue{1, 3}'],'Labels',{'Correct','Omature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 25])


%% 4. Get Z score normalized PSTHs for plotting colormap
        
% Get baseline subtracted PETHs
for n = 1:size(SD_OI.(SDperiod).psthBinsValue,2)
    OI_PSTH.PSTH.TS.Cor(n,:) = (mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:45)/binSize,1)-mean(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));;
    OI_PSTH.PSTH.TH.Cor(n,:) = (mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.wait_start(:,1:25)/binSize,1)-mean(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));
    OI_PSTH.PSTH.Resp.Cor(n,:) = (mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.resp(:,1:20)/binSize,1)-mean(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));

    OI_PSTH.PSTH.TS.Om(n,:) = (mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Om.trial_start(:,1:45)/binSize,1)-mean(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));
    OI_PSTH.PSTH.TH.Om(n,:) = (mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Om.wait_start(:,1:25)/binSize,1)-mean(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));
end

for n = 1:size(SD_All.(SDperiod).psthBinsValue,2)
    All_PSTH.PSTH.TS.Cor(n,:) = (mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:45)/binSize,1)-mean(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));
    All_PSTH.PSTH.TH.Cor(n,:) = (mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.wait_start(:,1:25)/binSize,1)-mean(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));
    All_PSTH.PSTH.Resp.Cor(n,:) = (mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.resp(:,1:20)/binSize,1)-mean(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));
    All_PSTH.PSTH.TS.Om(n,:) = (mean(SD_All.(SDperiod).psthBinsValue{1, n}.Om.trial_start(:,1:45)/binSize,1)-mean(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));
    All_PSTH.PSTH.TH.Om(n,:) = (mean(SD_All.(SDperiod).psthBinsValue{1, n}.Om.wait_start(:,1:25)/binSize,1)-mean(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));
end

time_axis_TS = linspace(-2,7,((2+7)/binSize));
time_axis_TH = linspace(-2,3,((2+3)/binSize));
time_axis_Resp = linspace(-2,2,((2+2)/binSize));
axisRange = [-6 6];

% Sort neurons in descending order on firing rate 2 sec before cue presentation per neuron catagory
[B,I_one] = sort(mean(OI_PSTH.PSTH.TS.Cor(find(ClustType_OI==1),26:35),2),'descend');
[B,I_zero] = sort(mean(OI_PSTH.PSTH.TS.Cor(find(ClustType_OI==0),26:35),2),'descend');
[B,I_neg] = sort(mean(OI_PSTH.PSTH.TS.Cor(find(ClustType_OI==-1),26:35),2),'descend');
I_oneIxs = find(ClustType_OI==1);
I_zeroIxs = find(ClustType_OI==0);
I_negIxs = find(ClustType_OI==-1);
I = [I_oneIxs(I_one) I_zeroIxs(I_zero) I_negIxs(I_neg)];

%% Figure 5E

figure
colormap('jet')
imagesc(B_OI(ib_OI(I),1:3))

%% Figure 5D

figure
subplot(1,3,1)
colormap('jet')
imagesc(time_axis_TS,1:(size(OI_PSTH.PSTH.TS.Cor,1)),OI_PSTH.PSTH.TS.Cor(I,1:45))
axis([-2-(0.5*binSize) 7+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.TS.Cor,1)+0.5])
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
%[B,I] = sort(mean(OI_PSTH.PSTH.cue_bef.Cor(:,1:10),2),'descend');
imagesc(time_axis_Resp,1:(size(OI_PSTH.PSTH.Resp.Cor,1)),OI_PSTH.PSTH.Resp.Cor(I,1:20))
axis([-2-(0.5*binSize) 2+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.Resp.Cor,1)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Resp')

%% Figure 5F

figure
subplot(1,3,1) 
hold on
SEM = nanstd(OI_PSTH.PSTH.TS.Cor(I_oneIxs(I_one),1:45),1)/sqrt(size(OI_PSTH.PSTH.TS.Cor(I_oneIxs(I_one),1:45),1)); 
shadedErrorBar(time_axis_TS,mean(OI_PSTH.PSTH.TS.Cor(I_oneIxs(I_one),1:45)),SEM,{'r'});
SEM = nanstd(OI_PSTH.PSTH.TS.Cor(I_zeroIxs(I_zero),1:45),1)/sqrt(size(OI_PSTH.PSTH.TS.Cor(I_zeroIxs(I_zero),1:45),1)); 
shadedErrorBar(time_axis_TS,mean(OI_PSTH.PSTH.TS.Cor(I_zeroIxs(I_zero),1:45)),SEM,{'k'});
SEM = nanstd(OI_PSTH.PSTH.TS.Cor(I_negIxs(I_neg),1:45),1)/sqrt(size(OI_PSTH.PSTH.TS.Cor(I_negIxs(I_neg),1:45),1)); 
shadedErrorBar(time_axis_TS,mean(OI_PSTH.PSTH.TS.Cor(I_negIxs(I_neg),1:45)),SEM,{'b'});
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
SEM = nanstd(OI_PSTH.PSTH.Resp.Cor(I_oneIxs(I_one),1:20),1)/sqrt(size(OI_PSTH.PSTH.Resp.Cor(I_oneIxs(I_one),1:20),1)); 
shadedErrorBar(time_axis_Resp,mean(OI_PSTH.PSTH.Resp.Cor(I_oneIxs(I_one),1:20)),SEM,{'r'});
SEM = nanstd(OI_PSTH.PSTH.Resp.Cor(I_zeroIxs(I_zero),1:20),1)/sqrt(size(OI_PSTH.PSTH.Resp.Cor(I_zeroIxs(I_zero),1:20),1)); 
shadedErrorBar(time_axis_Resp,mean(OI_PSTH.PSTH.Resp.Cor(I_zeroIxs(I_zero),1:20)),SEM,{'k'});
SEM = nanstd(OI_PSTH.PSTH.Resp.Cor(I_negIxs(I_neg),1:20),1)/sqrt(size(OI_PSTH.PSTH.Resp.Cor(I_negIxs(I_neg),1:20),1)); 
shadedErrorBar(time_axis_Resp,mean(OI_PSTH.PSTH.Resp.Cor(I_negIxs(I_neg),1:20)),SEM,{'b'});


I_oneIxs = find(ClustType_All(VAR_SDixOIMatched)==1);
I_zeroIxs = find(ClustType_All(VAR_SDixOIMatched)==0);
I_negIxs = find(ClustType_All(VAR_SDixOIMatched)==-1);
[B,I_one] = sort(mean(All_PSTH.PSTH.TS.Cor(VAR_SDixOIMatched(I_oneIxs),26:35),2),'descend');
[B,I_zero] = sort(mean(All_PSTH.PSTH.TS.Cor(VAR_SDixOIMatched(I_zeroIxs),26:35),2),'descend');
[B,I_neg] = sort(mean(All_PSTH.PSTH.TS.Cor(VAR_SDixOIMatched(I_negIxs),26:35),2),'descend');
I = [I_oneIxs(I_one) I_zeroIxs(I_zero) I_negIxs(I_neg)];

%% Figure S2B left

figure
colormap('jet')
imagesc(B_All(ib_All(VAR_SDixOIMatched(I)),1:3))

%% Figure S2A

figure
subplot(1,3,1)
colormap('jet')
imagesc(time_axis_TS,1:(size(All_PSTH.PSTH.TS.Cor(VAR_SDixOIMatched),2)),All_PSTH.PSTH.TS.Cor(VAR_SDixOIMatched(I),1:45))
axis([-2-(0.5*binSize) 7+(0.5*binSize) 0.5 size(All_PSTH.PSTH.TS.Cor(VAR_SDixOIMatched),2)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Trial start')

subplot(1,3,2)
colormap('jet')
imagesc(time_axis_TH,1:(size(All_PSTH.PSTH.TH.Cor(VAR_SDixOIMatched),2)),All_PSTH.PSTH.TH.Cor(VAR_SDixOIMatched(I),1:25))
axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(All_PSTH.PSTH.TH.Cor(VAR_SDixOIMatched),2)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Cue orientation')

subplot(1,3,3)
colormap('jet')
imagesc(time_axis_Resp,1:(size(All_PSTH.PSTH.Resp.Cor,1)),All_PSTH.PSTH.Resp.Cor(VAR_SDixOIMatched(I),1:20))
axis([-2-(0.5*binSize) 2+(0.5*binSize) 0.5 size(All_PSTH.PSTH.Resp.Cor,1)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Resp')
