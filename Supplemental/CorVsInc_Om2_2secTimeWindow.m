function CorVsInc_Om2_2secTimeWindow(All,neuronIndex, SD_All, SD_All, SD_OI, SD_OI, SDperiod)

% Structure
% 1. Get CorFr for OI and All per sync point. Calculate how much
% significant
% 2. Make scatterplots per sync point per up, down and unm
% 3. Get boxplots with signrank comp per sync point and up, down or unm
% 4. Get Z-score PSTHs for plotting colormap
% 5. Get Om cue times per SD type and plot
% 6. Get Cue orientation times
% 7. Get cueonse type distribution
binSize = 0.2;
SDperiod = 'All';
trialTypo = {'TS','TH','Cue'};
trialTypo2 = {'TS','TH','Cue'};
timeBins = {11:15, 48:57, 26:35};
%% 1. Get CorFr for OI and All per sync point

for i = 1:numel(neuronIndex.PyrIx_VAR_SD_OI)
    %ClustType_OI(i) = All{23,1}(neuronIndex.PyrIx_VAR_SD_OI(i),2);
    PSTH_Cor_OI{i} = [SD_OI.(SDperiod).psthBinsValue{1, i}.Cor.trial_start(:,1:35) SD_OI.(SDperiod).psthBinsValue{1, i}.Cor.wait_start(:,1:25)];
    PSTH_Om_OI{i} = [SD_OI.(SDperiod).psthBinsValue{1, i}.Om.trial_start(:,1:35) SD_OI.(SDperiod).psthBinsValue{1, i}.Om.wait_start(:,1:25)];
%     TresholdTime_Cor_OI(i) = mean(SD_OI.(SDperiod).TrialDist{1, i}.CorTreshold);
%     TresholdTime_Om_OI(i) = mean(SD_OI.(SDperiod).TrialDist{1, i}.OmTreshold);
%     TresholdTime.Cor_Raw_OI{i} = SD_OI.(SDperiod).TrialDist{1, i}.CorTreshold;
%     TresholdTime.Om_Raw_OI{i} = SD_OI.(SDperiod).TrialDist{1, i}.OmTreshold;
end

for i = 1:numel(neuronIndex.PyrIxVAR_SD)
    %ClustType_All(i) = All{23,1}(neuronIndex.PyrIxVAR_SD(i),2);
    PSTH_Cor_All{i} = [SD_All.(SDperiod).psthBinsValue{1, i}.Cor.trial_start(:,1:35) SD_All.(SDperiod).psthBinsValue{1, i}.Cor.wait_start(:,1:25)];
    PSTH_Om_All{i} = [SD_All.(SDperiod).psthBinsValue{1, i}.Om.trial_start(:,1:35) SD_All.(SDperiod).psthBinsValue{1, i}.Om.wait_start(:,1:25)];
%     TresholdTime_Cor_All(i) = mean(SD_All.(SDperiod).TrialDist{1, i}.CorTreshold);
%     TresholdTime_Om_All(i) = mean(SD_All.(SDperiod).TrialDist{1, i}.OmTreshold);
%     TresholdTime.Cor_Raw_All{i} = SD_All.(SDperiod).TrialDist{1, i}.CorTreshold;
%     TresholdTime.Om_Raw_All{i} = SD_All.(SDperiod).TrialDist{1, i}.OmTreshold;
end

%% 1. Get CorFr for OI and All per sync point, calculate significance

for t = 1:size(timeBins,2)
    for i = 1:numel(neuronIndex.PyrIx_VAR_SD_OI)

%             h1(i) = lillietest(mean(PSTH_Cor_OI{1, i}(:,timeBins{t}),2));
%             h2(i) = lillietest(mean(PSTH_Om_OI{1, i}(:,timeBins{t}),2)); 
% 
%             if h1(i) == 0 & h2(i) == 0
%                 [~, p] = ttest2(mean(PSTH_Cor_OI{1, i}(:,timeBins{t}),2)/0.2,mean(PSTH_Om_OI{1, i}(:,timeBins{t})/0.2,2),'Vartype','unequal');
%                 Sign.SD_OI.(trialTypo2{t})(i) = p;
%             else
%                 p = ranksum(mean(PSTH_Cor_OI{1, i}(:,timeBins{t}),2)/0.2,mean(PSTH_Om_OI{1, i}(:,timeBins{t})/0.2,2));
%                 Sign.SD_OI.(trialTypo2{t})(i) = p;
%             end
            p = ranksum(mean(PSTH_Cor_OI{1, i}(:,timeBins{t}),2)/binSize,mean(PSTH_Om_OI{1, i}(:,timeBins{t})/binSize,2));
            Sign.SD_OI.(trialTypo2{t})(i) = p;
            %CorFR.Raw.(trialTypo2{t}){i} = median(PSTH_Cor_OI{1, i}(:,timeBins{t})/0.2,2);
            %OmFR.Raw.(trialTypo2{t}){i} = median(PSTH_Om_OI{1, i}(:,timeBins{t})/0.2,2);
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
%             h1(i) = lillietest(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t}),2));
%             h2(i) = lillietest(mean(PSTH_Om_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t}),2)); 
%             if h1(i) == 0 & h2(i) == 0
%                 [~, p] = ttest2(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t}),2)/0.2,mean(PSTH_Om_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t}),2)/0.2,'Vartype','unequal');
%                 Sign.SD_All.(trialTypo2{t})(i) = p;
%             else
%                 p = ranksum(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t}),2)/0.2,mean(PSTH_Om_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t}),2)/0.2);
%                 Sign.SD_All.(trialTypo2{t})(i) = p;
%             end
            p = ranksum(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t}),2)/binSize,mean(PSTH_Om_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t}),2)/binSize);
            Sign.SD_All.(trialTypo2{t})(i) = p;
            CorFR.All.(trialTypo2{t})(i) = mean(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t})/binSize,2)-mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,1:10)/binSize,2));
            OmFR.All.(trialTypo2{t})(i) = mean(mean(PSTH_Om_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t})/binSize,2)-mean(PSTH_Om_All{1, VAR_SDixOIMatched(i)}(:,1:10)/binSize,2));
            CorFR_AUC.All.(trialTypo2{t})(i,:) = sum(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t})/0.2,1)-mean(mean(PSTH_Cor_All{1, VAR_SDixOIMatched(i)}(:,1:10)/0.2,1)));
            OmFR_AUC.All.(trialTypo2{t})(i,:) = sum(mean(PSTH_Om_All{1, VAR_SDixOIMatched(i)}(:,timeBins{t})/0.2,1)-mean(mean(PSTH_Om_All{1, VAR_SDixOIMatched(i)}(:,1:10)/0.2,1)));
    end
end
%Sign.SD_All =  mafdr([Sign.SD_All],'BHFDR',true);
%Sign.SD_OI =  mafdr([Sign.SD_OI],'BHFDR',true);


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
% 
% for i = 1:size(neuronIndex.PyrIx_VAR_SD_OI,2)
%     if ib_OI(i) == 1 || ib_OI(i) == 5 || ib_OI(i) == 12
%         ClustType_OI(i) = -1;
%     elseif ib_OI(i) == 8
%         ClustType_OI(i) = 0;
%     elseif ib_OI(i) == 15 || ib_OI(i) == 11 || ib_OI(i) == 4
%         ClustType_OI(i) = 1;
%     else
%         ClustType_OI(i) = 0;
%     end
% end

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

% for i = 1:size(neuronIndex.PyrIxVAR_SD,2)
%     if ib_All(i) == 1 || ib_All(i) == 8 || ib_All(i) == 16
%         ClustType_All(i) = -1;
%     elseif ib_All(i) == 12
%         ClustType_All(i) = 0;
%     elseif ib_All(i) == 23 || ib_All(i) == 15 || ib_All(i) == 7
%         ClustType_All(i) = 1;
%     else
%         ClustType_All(i) = 0;
%     end
% end

SD_OI_clust.TS = ClustType_OI';
SD_OI_clust.TH = ClustType_OI';
SD_OI_clust.cue_bef = ClustType_OI';

SD_All_clust.TS = ClustType_All';
SD_All_clust.TH = ClustType_All';
SD_All_clust.cue_bef = ClustType_All';

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
%Chi_cueTypesOIvsOther.pval = mafdr(Chi_cueTypesOIvsOther.pval,'BHFDR',true);


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
            %subplot(4,3,pos(n))
            hold on
            %plot([CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))); OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))],'k')
            %plot([CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)& Sign.SD_OI.(trialTypo2{t}) < 0.05)); OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl) & Sign.SD_OI.(trialTypo2{t}) < 0.05))],'r')
            %boxplot([CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))' OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))'],'Labels',{'Correct','Omature'})
            diffCorOm.OI.(trialTypo2{t})(Cl) = median(((OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))'-CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))')./abs((CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))')))*100);
            %ylabel('Firing rate (Hz)')
            %diffCorOm.OI2.(trialTypo2{t}){Cl} = ((OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))'-CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))')./(CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))'))*100;

            try
            CorVsOmSign.OI.(trialTypo2{t})(Cl) = signrank(CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))));
            %CorVsOmSign.OI.(trialTypo2{t})(Cl) = signrank(CorFR_AUC.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),OmFR_AUC.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))));

            %CorVsOmSign.OIVal.(trialTypo2{t}){Cl} = abs(CorFR.OI.(trialTypo2{t})(find(Sign.OI.(trialTypo2{t}) < 0.05))-OmFR.OI.(trialTypo2{t})(find(Sign.OI.(trialTypo2{t}) < 0.05)));
            catch
                CorVsOmSign.OI.(trialTypo2{t})(Cl) = nan;
            end
        end
        subplot(3,1,pos(n))
        bar(diffCorOm.OI.(trialTypo2{t}))
        ylim([-25 25])
        n = n+1;
        %CorVsOmSign.OI.(trialTypo2{t}) = mafdr([CorVsOmSign.OI.(trialTypo2{t})],'BHFDR',true);
        CorVsOmSign_Combined.OI.(trialTypo2{t}) = signrank(abs(CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))),abs(OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))));
        CorVsOmSign_Combined2.OI.(trialTypo2{t}){2} = abs(CorFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2))));
        CorVsOmSign_Combined2.OI.(trialTypo2{t}){3} = abs(OmFR.OI.(trialTypo2{t})(find(SD_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2))));
end
%CorVsOmSign_Combined.OI_combined = mafdr([CorVsOmSign_Combined.OI.TH CorVsOmSign_Combined.OI.Cue],'BHFDR',true);
CorVsOmSign_Combined.OI_combined = [CorVsOmSign_Combined.OI.TH CorVsOmSign_Combined.OI.Cue];
CorVsOmSign_Combined.OI.TH = CorVsOmSign_Combined.OI_combined(1);
CorVsOmSign_Combined.OI.Cue = CorVsOmSign_Combined.OI_combined(2);
figure
subplot(1,3,1)
hold on
plot([abs(CorFR.OI.TS(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1)))); abs(OmFR.OI.TS(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1))))],'r');
plot([abs(CorFR.OI.TS(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3)))); abs(OmFR.OI.TS(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3))))],'b');
boxplot([CorVsOmSign_Combined2.OI.TS{1, 2}' CorVsOmSign_Combined2.OI.TS{1, 3}'],'Labels',{'Correct','Omature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 20])
subplot(1,3,2)
hold on
plot([abs(CorFR.OI.TH(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1)))); abs(OmFR.OI.TH(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1))))],'r');
plot([abs(CorFR.OI.TH(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3)))); abs(OmFR.OI.TH(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3))))],'b');
boxplot([CorVsOmSign_Combined2.OI.TH{1, 2}' CorVsOmSign_Combined2.OI.TH{1, 3}'],'Labels',{'Correct','Omature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 20])
subplot(1,3,3)
hold on
plot([abs(CorFR.OI.Cue(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1)))); abs(OmFR.OI.Cue(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(1))))],'r');
plot([abs(CorFR.OI.Cue(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3)))); abs(OmFR.OI.Cue(find(SD_OI_clust.(trialTypo2{t})'==ClustUpMidDown(3))))],'b');
boxplot([CorVsOmSign_Combined2.OI.Cue{1, 2}' CorVsOmSign_Combined2.OI.Cue{1, 3}'],'Labels',{'Correct','Omature'})
ylabel('Absolute delta FR (Hz)')
ylim([0 20])

figure
n = 1;
pos = [1 2 3];
for t = [1 2 3]
        for Cl = 1:numel(ClustUpMidDown)
            %subplot(4,3,pos(n))
            %hold on
            %plot([CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))); OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))],'k')
            %plot([CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)& Sign.SD_All.(trialTypo2{t}) < 0.05)); OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl) & Sign.SD_All.(trialTypo2{t}) < 0.05))],'r')
            %boxplot([CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))' OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))'],'Labels',{'Correct','Omature'})
            %ylabel('Firing rate (Hz)')
            diffCorOm.All.(trialTypo2{t})(Cl) = median(((OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))'-CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))')./abs((CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))')))*100);
            %diffCorOm.All2.(trialTypo2{t}){Cl} = ((OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))'-CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))')./(CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))'))*100;

            %n = n+1;
            h1(i) = lillietest(CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))));
            h2(i) = lillietest(OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl)))); 
            if h1(i) == 0 & h2(i) == 0
                try
                [h CorVsOmSign.All.(trialTypo2{t})(Cl)] = ttest(CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))),OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))),'Vartype','unequal');
                %CorVsOmSign.AllVal.(trialTypo2{t}){Cl} = abs(CorFR.All.(trialTypo2{t})(find(Sign.SD_All.(trialTypo2{t}) < 0.05))-OmFR.All.(trialTypo2{t})(find(Sign.SD_All.(trialTypo2{t}) < 0.05)));

                catch
                    CorVsOmSign.All.(trialTypo2{t})(Cl) = nan;
                end
            else
                try
                CorVsOmSign.All.(trialTypo2{t})(Cl) = signrank(CorFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))),OmFR.All.(trialTypo2{t})(find(SD_All_clust.(trialTypo2{t})(VAR_SDixOIMatched)'==ClustUpMidDown(Cl))));
                %CorVsOmSign.AllVal.(trialTypo2{t}){Cl} = abs(CorFR.All.(trialTypo2{t})(find(Sign.SD_All.(trialTypo2{t}) < 0.05))-OmFR.All.(trialTypo2{t})(find(Sign.SD_All.(trialTypo2{t}) < 0.05)));

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
%CorVsOmSign_Combined.All_combined = mafdr([CorVsOmSign_Combined.All.TH CorVsOmSign_Combined.All.Cue],'BHFDR',true);
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
%timeBins = {11:15, 11:20, 1:10, 12:20, 1:10};

        
% Get baseline subtracted PETHs
for n = 1:size(SD_OI.(SDperiod).psthBinsValue,2)
    OI_PSTH.PSTH.TS.Cor(n,:) = (mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:45)/binSize,1)-mean(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));% /max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,11:15),1)));
    OI_PSTH.PSTH.TH.Cor(n,:) = (mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.wait_start(:,1:25)/binSize,1)-mean(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));%/max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.wait_start(:,6:25),1)));
    OI_PSTH.PSTH.Resp.Cor(n,:) = (mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.resp(:,1:20)/binSize,1)-mean(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));%/max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.wait_start(:,6:25),1)));

    %OI_PSTH.PSTH.cue_bef.Cor(n,:) = (mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.cue(:,1:30)/binSize,1)-mean(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));%/max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.cue(:,1:10),1)));
    OI_PSTH.PSTH.TS.Om(n,:) = (mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Om.trial_start(:,1:45)/binSize,1)-mean(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Om.trial_start(:,1:10)/binSize,1));%/max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,11:15),1)));
    OI_PSTH.PSTH.TH.Om(n,:) = (mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Om.wait_start(:,1:25)/binSize,1)-mean(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Om.trial_start(:,1:10)/binSize,1));%/max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.wait_start(:,6:25),1)));
    %OI_PSTH.PSTH.cue_bef.Om(n,:) = (mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Om.cue(:,1:30)/binSize,1)-mean(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Om.trial_start(:,1:10)/binSize,1));%/max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Om.cue(:,1:10),1)));
end

for n = 1:size(SD_All.(SDperiod).psthBinsValue,2)
    All_PSTH.PSTH.TS.Cor(n,:) = (mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:45)/binSize,1)-mean(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));% /max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,11:15),1)));
    All_PSTH.PSTH.TH.Cor(n,:) = (mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.wait_start(:,1:25)/binSize,1)-mean(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));%/max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.wait_start(:,6:25),1)));
    All_PSTH.PSTH.Resp.Cor(n,:) = (mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.resp(:,1:20)/binSize,1)-mean(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));%/max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.wait_start(:,6:25),1)));

    %All_PSTH.PSTH.cue_bef.Cor(n,:) = (mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.cue(:,1:30)/binSize,1)-mean(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1));%/max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.cue(:,1:10),1)));
    All_PSTH.PSTH.TS.Om(n,:) = (mean(SD_All.(SDperiod).psthBinsValue{1, n}.Om.trial_start(:,1:45)/binSize,1)-mean(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Om.trial_start(:,1:10)/binSize,1));%/max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,11:15),1)));
    All_PSTH.PSTH.TH.Om(n,:) = (mean(SD_All.(SDperiod).psthBinsValue{1, n}.Om.wait_start(:,1:25)/binSize,1)-mean(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Om.trial_start(:,1:10)/binSize,1));%/max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Cor.wait_start(:,6:25),1)));
    %All_PSTH.PSTH.cue_bef.Om(n,:) = (mean(SD_All.(SDperiod).psthBinsValue{1, n}.Om.cue(:,1:30)/binSize,1)-mean(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Cor.trial_start(:,1:10)/binSize,1)));%/std(mean(SD_All.(SDperiod).psthBinsValue{1, n}.Om.trial_start(:,1:10)/binSize,1));%/max(abs(mean(SD_OI.(SDperiod).psthBinsValue{1, n}.Om.cue(:,1:10),1)));
end

time_axis_TS = linspace(-2,7,((2+7)/binSize));
time_axis_TH = linspace(-2,3,((2+3)/binSize));
time_axis_Resp = linspace(-2,2,((2+2)/binSize));
%time_axis_cue_bef = linspace(-2,4,((2+4)/binSize));
axisRange = [-6 6];

[B,I_one] = sort(mean(OI_PSTH.PSTH.TS.Cor(find(ClustType_OI==1),26:35),2),'descend');
[B,I_zero] = sort(mean(OI_PSTH.PSTH.TS.Cor(find(ClustType_OI==0),26:35),2),'descend');
[B,I_neg] = sort(mean(OI_PSTH.PSTH.TS.Cor(find(ClustType_OI==-1),26:35),2),'descend');
I_oneIxs = find(ClustType_OI==1);
I_zeroIxs = find(ClustType_OI==0);
I_negIxs = find(ClustType_OI==-1);
I = [I_oneIxs(I_one) I_zeroIxs(I_zero) I_negIxs(I_neg)];

figure
colormap('jet')
imagesc(B_OI(ib_OI(I),1:3))

% [B,I] = sort(mean(OI_PSTH.PSTH.TH.Cor(:,13:17),2),'descend');
figure
subplot(1,3,1)
colormap('jet')
%[B,I] = sort(mean(OI_PSTH.PSTH.TS.Cor(:,11:15),2),'descend');
imagesc(time_axis_TS,1:(size(OI_PSTH.PSTH.TS.Cor,1)),OI_PSTH.PSTH.TS.Cor(I,1:45))
axis([-2-(0.5*binSize) 7+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.TS.Cor,1)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Trial start')

% subplot(2,3,4)
% colormap('jet')
% imagesc(time_axis_TS,1:(size(OI_PSTH.PSTH.TS.Om,1)),OI_PSTH.PSTH.TS.Om(I,1:45))
% axis([-2-(0.5*binSize) 7+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.TS.Om,1)+0.5])
% ylabel('Neuron')
% caxis(axisRange)
% title('Trial start')

subplot(1,3,2)
colormap('jet')
%[B,I] = sort(mean(OI_PSTH.PSTH.TH.Cor(:,11:15),2),'descend');
imagesc(time_axis_TH,1:(size(OI_PSTH.PSTH.TH.Cor,1)),OI_PSTH.PSTH.TH.Cor(I,1:25))
axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.TH.Cor,1)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Cue orientation')

% subplot(2,3,5)
% colormap('jet')
% imagesc(time_axis_TH,1:(size(OI_PSTH.PSTH.TH.Om,1)),OI_PSTH.PSTH.TH.Om(I,1:25))
% axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.TH.Om,1)+0.5])
% ylabel('Neuron')
% caxis(axisRange)
% title('Cue orientation')

subplot(1,3,3)
colormap('jet')
%[B,I] = sort(mean(OI_PSTH.PSTH.cue_bef.Cor(:,1:10),2),'descend');
imagesc(time_axis_Resp,1:(size(OI_PSTH.PSTH.Resp.Cor,1)),OI_PSTH.PSTH.Resp.Cor(I,1:20))
axis([-2-(0.5*binSize) 2+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.Resp.Cor,1)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Resp')

% subplot(2,3,3)
% colormap('jet')
% %[B,I] = sort(mean(OI_PSTH.PSTH.cue_bef.Cor(:,1:10),2),'descend');
% imagesc(time_axis_cue_bef,1:(size(OI_PSTH.PSTH.cue_bef.Cor,1)),OI_PSTH.PSTH.Cue.Cor(I,1:30))
% axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.Cue.Cor,1)+0.5])
% ylabel('Neuron')
% caxis(axisRange)
% title('cueonse')
% 
% subplot(2,3,6)
% colormap('jet')
% imagesc(time_axis_cue_bef,1:(size(OI_PSTH.PSTH.cue_bef.Om,1)),OI_PSTH.PSTH.Cue.Om(I,1:30))
% axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.Cue.Om,1)+0.5])
% ylabel('Neuron')
% caxis(axisRange)
% title('cueonse')

% Plot Average PSTH per group
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

figure
colormap('jet')
imagesc(B_All(ib_All(VAR_SDixOIMatched(I)),1:3))

figure
subplot(1,3,1)
colormap('jet')
%[B,I] = sort(mean(All_PSTH.PSTH.TS.Cor(VAR_SDixOIMatched,11:15),2),'descend');
imagesc(time_axis_TS,1:(size(All_PSTH.PSTH.TS.Cor(VAR_SDixOIMatched),2)),All_PSTH.PSTH.TS.Cor(VAR_SDixOIMatched(I),1:45))
axis([-2-(0.5*binSize) 7+(0.5*binSize) 0.5 size(All_PSTH.PSTH.TS.Cor(VAR_SDixOIMatched),2)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Trial start')

% subplot(2,3,4)
% colormap('jet')
% imagesc(time_axis_TS,1:(size(All_PSTH.PSTH.TS.Om(VAR_SDixOIMatched),2)),All_PSTH.PSTH.TS.Om(VAR_SDixOIMatched(I),1:45))
% axis([-2-(0.5*binSize) 7+(0.5*binSize) 0.5 size(All_PSTH.PSTH.TS.Om(VAR_SDixOIMatched),2)+0.5])
% ylabel('Neuron')
% caxis(axisRange)
% title('Trial start')

subplot(1,3,2)
colormap('jet')
%[B,I] = sort(mean(All_PSTH.PSTH.TH.Cor(VAR_SDixOIMatched,11:15),2),'descend');
imagesc(time_axis_TH,1:(size(All_PSTH.PSTH.TH.Cor(VAR_SDixOIMatched),2)),All_PSTH.PSTH.TH.Cor(VAR_SDixOIMatched(I),1:25))
axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(All_PSTH.PSTH.TH.Cor(VAR_SDixOIMatched),2)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Cue orientation')

% subplot(2,3,5)
% colormap('jet')
% imagesc(time_axis_TH,1:(size(All_PSTH.PSTH.TH.Om(VAR_SDixOIMatched),2)),All_PSTH.PSTH.TH.Om(VAR_SDixOIMatched(I),1:25))
% axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(All_PSTH.PSTH.TH.Om(VAR_SDixOIMatched),2)+0.5])
% ylabel('Neuron')
% caxis(axisRange)
% title('Cue orientation')

subplot(1,3,3)
colormap('jet')
%[B,I] = sort(mean(OI_PSTH.PSTH.cue_bef.Cor(:,1:10),2),'descend');
imagesc(time_axis_Resp,1:(size(All_PSTH.PSTH.Resp.Cor,1)),All_PSTH.PSTH.Resp.Cor(VAR_SDixOIMatched(I),1:20))
axis([-2-(0.5*binSize) 2+(0.5*binSize) 0.5 size(All_PSTH.PSTH.Resp.Cor,1)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Resp')

% subplot(2,3,3)
% colormap('jet')
% %[B,I] = sort(mean(All_PSTH.PSTH.cue_bef.Cor(VAR_SDixOIMatched,1:10),2),'descend');
% imagesc(time_axis_cue_bef,1:(size(All_PSTH.PSTH.cue_bef.Cor(VAR_SDixOIMatched),2)),All_PSTH.PSTH.cue_bef.Cor(VAR_SDixOIMatched(I),1:30))
% axis([-2-(0.5*binSize) 4+(0.5*binSize) 0.5 size(All_PSTH.PSTH.cue_bef.Cor(VAR_SDixOIMatched),2)+0.5])
% ylabel('Neuron')
% caxis(axisRange)
% title('cueonse')
% 
% subplot(2,3,6)
% colormap('jet')
% imagesc(time_axis_cue_bef,1:(size(All_PSTH.PSTH.cue_bef.Om(VAR_SDixOIMatched),2)),All_PSTH.PSTH.cue_bef.Om(VAR_SDixOIMatched(I),1:30))
% axis([-2-(0.5*binSize) 4+(0.5*binSize) 0.5 size(All_PSTH.PSTH.cue_bef.Om(VAR_SDixOIMatched),2)+0.5])
% ylabel('Neuron')
% caxis(axisRange)
% title('cueonse')

%% Make plot that compares proportion of Cor vs Om modualted units per sync point for OI and All

PropTSMod_OI = numel(find(Sign.SD_OI.TS<0.05))/(numel(Sign.SD_OI.TS));
PropTHMod_OI  = numel(find(Sign.SD_OI.TH<0.05))/(numel(Sign.SD_OI.TH));
Propcue_befMod_OI  = numel(find(Sign.SD_OI.Cue<0.05))/(numel(Sign.SD_OI.Cue));
PropTSMod_All = numel(find(Sign.SD_All.TS<0.05))/(numel(Sign.SD_All.TS));
PropTHMod_All  = numel(find(Sign.SD_All.TH<0.05))/(numel(Sign.SD_All.TH));
Propcue_befMod_All  = numel(find(Sign.SD_All.Cue<0.05))/(numel(Sign.SD_All.Cue));
figure
bar([PropTSMod_OI PropTSMod_All;PropTHMod_OI PropTHMod_All; Propcue_befMod_OI Propcue_befMod_All])
title('Proportion of Cor-Om modulated neurons')
legend({'Identified','other'})

% ClustUpMidDown = [1 0 -1];
% tot_All_OI = numel(SD_OI_clust.TS);
% tot_All_All = numel(VAR_SDixOIMatched);
% for i = 1:numel(ClustUpMidDown)
%     x = [zeros(numel(SD_OI_clust.TS),1); ones(numel(VAR_SDixOIMatched),1)]';
%     y = [zeros(numel(find(SD_OI_clust.TS == ClustUpMidDown(i))),1); ones(tot_All_OI-numel(find(SD_OI_clust.TS == ClustUpMidDown(i))),1); zeros(numel(find(SD_All_clust.TS(VAR_SDixOIMatched) == ClustUpMidDown(i))),1); ones(tot_All_All-numel(find(SD_All_clust.TS(VAR_SDixOIMatched) == ClustUpMidDown(i))),1)]';
%     [tbl_cueTypesOIvsOther{i},chi2stat_cueTypesOIvsOther{i},Chi_cueTypesOIvsOther.pval(i)] = crosstab(x,y);
%     if numel(find(SD_OI_clust.TS == ClustUpMidDown(i))) < 5 || numel(find(SD_All_clust.TS(VAR_SDixOIMatched) == ClustUpMidDown(i))) < 5
%         tbl_cueTypesOIvsOther{i} = table([numel(find(SD_OI_clust.TS == ClustUpMidDown(i)));tot_All_OI-numel(find(SD_OI_clust.TS == ClustUpMidDown(i)))],[numel(find(SD_All_clust.TS(VAR_SDixOIMatched) == ClustUpMidDown(i)));tot_All_All-numel(find(SD_All_clust.TS(VAR_SDixOIMatched) == ClustUpMidDown(i)))],'VariableNames',{'Flu','NoFlu'},'RowNames',{'NoShot','Shot'});
%         [~,Chi_cueTypesOIvsOther.pval(i),chi2stat_cueTypesOIvsOther{i}] = fishertest(tbl_cueTypesOIvsOther{i});
%     end
% end
% Chi_cueTypesOIvsOther.pval = mafdr(Chi_cueTypesOIvsOther.pval,'BHFDR',true);

%% Plot Om cue times and Cue orientation times
allFiles = dir;

n = 1;
VAR_SDix = [];
VAR_SDix = [];
for file = 1:size(allFiles,1)
    if ~isempty(strfind(allFiles(file).name, 'txt'))
        MED_PC{n,1} = allFiles(file).name;
        n = n+1;
    end
end

for txtFile = 1:numel(MED_PC)
    for i = 1:size(All{5, 1},1)
        if All{5, 1}{i,2} == MED_PC{txtFile,1} & contains(All{5, 1}{i,1}, 'SD')
            VAR_SDix(end+1) = txtFile;
        elseif All{5, 1}{i,2} == MED_PC{txtFile,1} & contains(All{5, 1}{i,1}, 'SD')
            VAR_SDix(end+1) = txtFile;
        end
    end
end

VAR_SDix = unique(VAR_SDix);


% Get performance per session
n = 1;
for sessionSD = VAR_SDix
    [SD.Correct{n,1}, SD.Incorrect{n,1}, SD.Omission{n,1}, OmatureSD(n,1), SD.Perserverative{n,1}, nrTrialsStageSD(n,1), Perf_SD{n,1}]  = CC_ephys_var_ITI_SD_TTL_Master(10,1,MED_PC{sessionSD,1});
    n = n+1;
end

Lat.SD50 = [];
Lat.SD75 = [];
Lat.SD125 = [];



n = 1;
for i = 1:size(OmatureSD,1)
    Lat.SD125(n,:) = histcounts((OmatureSD(i).Latency),0:0.5:5)/numel(OmatureSD(i).Latency);
    n = n+1;
end
for i = 1:size(Lat.SD125,1)
    if Lat.SD50(i,11) > 0
        Lat.SD50(i,10) = Lat.SD125(i,10)+Lat.SD125(i,11);
        Lat.SD50(i,11) = 0;
    elseif Lat.SD75(i,16) > 0
        Lat.SD75(i,15) = Lat.SD125(i,15)+Lat.SD125(i,16);
        Lat.SD75(i,16) = 0;
    end
end

n = 2;
TH_Time_Om{1} = SD_OI.(SDperiod).TrialDist{1, 1}.OmTreshold;
TH_Time_cor{1} = SD_OI.(SDperiod).TrialDist{1, 1}.CorTreshold;
for neuron = 2:size(SD_OI.(SDperiod).TrialDist,2)
    if numel(SD_OI.(SDperiod).TrialDist{1, neuron}.OmTreshold)~=numel(TH_Time_Om{n-1})
        TH_Time_Om{n} = SD_OI.(SDperiod).TrialDist{1, neuron}.OmTreshold;
        TH_Time_cor{n} = SD_OI.(SDperiod).TrialDist{1, neuron}.CorTreshold;
        n = n+1;
    end
end

n = 1;
for i = 1:size(TH_Time_Om,2)
    Lat_TH_Om(n,:) = histcounts((TH_Time_Om{i}),0:0.5:5)/numel(TH_Time_Om{i});
    Lat_TH_cor(n,:) = histcounts((TH_Time_cor{i}),0:0.5:5)/numel(TH_Time_cor{i});
    n = n+1;
end

n = 1;
for i = 1:size(TH_Time_Om,2)
    Lat_TH_Om2(n) = median(TH_Time_Om{i});
    Lat_TH_cor2(n) = median(TH_Time_cor{i});
    n = n+1;
end

x_SD125 = 0.5:0.5:5;

figure
plot(x_SD125,[mean(Lat.SD125,1)'])
hold on
plot(x_SD125,[mean(Lat_TH_Om,1)'])
plot(x_SD125,[mean(Lat_TH_cor,1)'])
title('Omature cueonse time Var. SD')
xlabel('Time (sec)')
ylabel('cueonses/session (mean)')
legend({'SD 5.0'})

figure
errorbar(x_SD125,[mean(Lat_TH_cor,1)'],std(Lat_TH_cor,1))
title('Omature cueonse time Var. SD')
xlabel('Time (sec)')
ylabel('cueonses/session (mean)')
legend({'SD 5.0'})

figure
hold on
errorbar(x_SD125,[mean(Lat_TH_cor,1)'],std(Lat_TH_cor,1))
errorbar(x_SD125,[mean(Lat_TH_Om,1)'],std(Lat_TH_Om,1))
title('Omature cueonse time Var. SD')
xlabel('Time (sec)')
ylabel('cueonses/session (mean)')
%legend({'SD 5.0'})

Lat_TH_cor2(Lat_TH_cor2<0.1) = [] ;
Lat_TH_Om2(Lat_TH_Om2<0.1) = [];
figure
hold on
boxplot([Lat_TH_cor2' Lat_TH_Om2'],'Labels',{'Correct','Omission'})
scatter(ones(size(Lat_TH_cor2,2),1).*(1+(rand(size(Lat_TH_cor2,2),1)-0.5)/5),Lat_TH_cor2,20,'k','filled')
scatter(ones(size(Lat_TH_Om2,2),1).*(2+(rand(size(Lat_TH_Om2,2),1)-0.5)/5),Lat_TH_Om2,20,'k','filled')
ylabel('Start trial - Cue orientation (sec)')
ylim([0 5])

%% Junk

for t = [1 2 5]
    tot_All_OI = 
    tot_All_All
    for Cl = 1:numel(ClustUpMidDown)
        x = [zeros(tot_All_OI,1); ones(tot_All_All,1)]';
        y = [zeros(numel(find([SD_OI_clust] == ClustUpMidDown(i))),1); ones(tot_All_OI-numel(find([SD_OI_clust] == ClustUpMidDown(i))),1); zeros(numel(find([SD_All_clust] == ClustUpMidDown(i))),1); ones(tot_All_All-numel(find([SD_All_clust] == ClustUpMidDown(i))),1)]';
        [tbl_All{i},chi2stat_All{i},Chi_All.pval(i)] = crosstab(x,y);
        if numel(find([SD_OI_clust] == ClustUpMidDown(i))) < 5 || numel(find([SD_All_clust] == ClustUpMidDown(i))) < 5
            tbl_All{i} = table([numel(find([SD_OI_clust] == ClustUpMidDown(i)));tot_All_OI-numel(find([SD_OI_clust] == ClustUpMidDown(i)))],[numel(find([SD_All_clust] == ClustUpMidDown(i)));tot_All_All-numel(find([SD_All_clust] == ClustUpMidDown(i)))],'VariableNames',{'Flu','NoFlu'},'RowNames',{'NoShot','Shot'});
            [~,Chi_All.pval(i),chi2stat_All{i}] = fishertest(tbl_All{i});
        end
    end
Chi_All.pval = mafdr(Chi_All.pval,'BHFDR',true);



%% Omature modulation index
for i = 1:size(CorFR.OI.cue_bef,2)
    AMI_OI.cue_bef(i) = (CorFR.OI.cue_bef(i)-OmFR.OI.cue_bef(i))/(CorFR.OI.cue_bef(i)+OmFR.OI.cue_bef(i));
end
for i = 1:size(CorFR.All.cue_bef,2)
    AMI_All.cue_bef(i) = (CorFR.All.cue_bef(i)-OmFR.All.cue_bef(i))/(CorFR.All.cue_bef(i)+OmFR.All.cue_bef(i));
end

for n = 1:size(SD_OI.All.PSTH,2)
    All_PSTH.TS.Cor(n,:) = (SD_OI.(SDperiod).PSTH{1, n}.Cor.trial_start-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:10)))/max(abs(SD_OI.All.PSTH{1, n}.Cor.trial_start-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:10))));%,'gaussian',3)/SD_OI.All.ID{1, n}.Cor.TrialStart.Peak_real;
    All_PSTH.TH.Cor(n,:) = (SD_OI.(SDperiod).PSTH{1, n}.Cor.wait_start-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:10)))/max(abs(SD_OI.All.PSTH{1, n}.Cor.wait_start-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:10))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Treshold.Peak_real;
    All_PSTH.cue.Cor(n,:) = (SD_OI.(SDperiod).PSTH{1, n}.Cor.cue-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:10)))/max(abs(SD_OI.All.PSTH{1, n}.Cor.cue-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:10))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.cueonse.Peak_real;
    All_PSTH.TS.Om(n,:) = (SD_OI.(SDperiod).PSTH{1, n}.Om.trial_start-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:10)))/max(abs(SD_OI.All.PSTH{1, n}.Cor.trial_start-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:10))));%,'gaussian',3)/SD_OI.All.ID{1, n}.Cor.TrialStart.Peak_real;
    All_PSTH.TH.Om(n,:) = (SD_OI.(SDperiod).PSTH{1, n}.Om.wait_start-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:10)))/max(abs(SD_OI.All.PSTH{1, n}.Cor.wait_start-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:10))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Treshold.Peak_real;
    All_PSTH.cue.Om(n,:) = (SD_OI.(SDperiod).PSTH{1, n}.Om.cue-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:10)))/max(abs(SD_OI.All.PSTH{1, n}.Cor.cue-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:10))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.cueonse.Peak_real;
end

for Cl = 1:numel(ClustUpMidDown)
subplot(2,3,Cl+3)
time_x = -2:binSize:12.8;
plot(time_x,[mean(All_PSTH.TS.Cor(find(Sign.SD_OI(1,:)<0.05&ClustType_OI(1,:)==ClustUpMidDown(Cl)),1:75),1)],'k')
hold on
plot(time_x,[mean(All_PSTH.TS.Om(find(Sign.SD_OI(1,:)<0.05&ClustType_OI(1,:)==ClustUpMidDown(Cl)),1:75),1)],'b')
legend({'Cor','Om'})
xlabel('Time cue ori')
ylabel('FR (cor max normalized)')
end