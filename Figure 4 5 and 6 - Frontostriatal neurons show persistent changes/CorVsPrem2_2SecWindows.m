function CorVsPrem2_2SecWindows(All,neuronIndex, SD_All, ITI_All, SD_OI, ITI_OI, ITIperiod)

% Structure
% 1. Get CorFr for OI and All per sync point. Calculate how much
% significant
% 2. Make scatterplots per sync point per up, down and unm
% 3. Get boxplots with signrank comp per sync point and up, down or unm
% 4. Get Z-score PSTHs for plotting colormap
% 5. Get Prem resp times per ITI type and plot
% 6. Get Cue orientation times
% 7. Get Response type distribution
binSize = 0.2;
ITIperiod = 'ITI125';
trialTypo = {'TS','TH','Resp'};
trialTypo2 = {'TS','TH','Resp_bef'};
%timeBins = {11:15, 38:42, 51:55};
timeBins = {11:15, 38:47, 51:60};
%% 1. Get CorFr for OI and All per sync point

for i = 1:numel(neuronIndex.PyrIx_VAR_ITI_OI)
    %ClustType_OI(i) = All{23,1}(neuronIndex.PyrIx_VAR_ITI_OI(i),2);
    PSTH_Cor_OI{i} = [ITI_OI.(ITIperiod).psthBinsValue{1, i}.Cor.trial_start(:,1:25) ITI_OI.(ITIperiod).psthBinsValue{1, i}.Cor.wait_start(:,1:25) ITI_OI.(ITIperiod).psthBinsValue{1, i}.Cor.resp(:,1:20)];
    PSTH_Prem_OI{i} = [ITI_OI.(ITIperiod).psthBinsValue{1, i}.Prem.trial_start(:,1:25) ITI_OI.(ITIperiod).psthBinsValue{1, i}.Prem.wait_start(:,1:25) ITI_OI.(ITIperiod).psthBinsValue{1, i}.Prem.resp(:,1:20)];
    TresholdTime_Cor_OI(i) = mean(ITI_OI.(ITIperiod).TrialDist{1, i}.CorTreshold);
    TresholdTime_Prem_OI(i) = mean(ITI_OI.(ITIperiod).TrialDist{1, i}.PremTreshold);
    TresholdTime.Cor_Raw_OI{i} = ITI_OI.(ITIperiod).TrialDist{1, i}.CorTreshold;
    TresholdTime.Prem_Raw_OI{i} = ITI_OI.(ITIperiod).TrialDist{1, i}.PremTreshold;
end

for i = 1:numel(neuronIndex.PyrIxVAR_ITI)
    %ClustType_All(i) = All{23,1}(neuronIndex.PyrIxVAR_ITI(i),2);
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

%             h1(i) = lillietest(mean(PSTH_Cor_OI{1, i}(:,timeBins{t}),2));
%             h2(i) = lillietest(mean(PSTH_Prem_OI{1, i}(:,timeBins{t}),2)); 
% 
%             if h1(i) == 0 & h2(i) == 0
%                 [~, p] = ttest2(mean(PSTH_Cor_OI{1, i}(:,timeBins{t}),2)/binSize,mean(PSTH_Prem_OI{1, i}(:,timeBins{t})/binSize,2),'Vartype','unequal');
%                 Sign.ITI_OI.(trialTypo2{t})(i) = p;
%             else
                p = ranksum(mean(PSTH_Cor_OI{1, i}(:,timeBins{t}),2)/binSize,mean(PSTH_Prem_OI{1, i}(:,timeBins{t})/binSize,2));
                Sign.ITI_OI.(trialTypo2{t})(i) = p;
%             end

            %CorFR.Raw.(trialTypo2{t}){i} = median(PSTH_Cor_OI{1, i}(:,timeBins{t})/binSize,2);
            %PremFR.Raw.(trialTypo2{t}){i} = median(PSTH_Prem_OI{1, i}(:,timeBins{t})/binSize,2);
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
%             h1(i) = lillietest(mean(PSTH_Cor_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t}),2));
%             h2(i) = lillietest(mean(PSTH_Prem_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t}),2)); 
%             if h1(i) == 0 & h2(i) == 0
%                 [~, p] = ttest2(mean(PSTH_Cor_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t}),2)/binSize,mean(PSTH_Prem_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t}),2)/binSize,'Vartype','unequal');
%                 Sign.ITI_All.(trialTypo2{t})(i) = p;
%             else
                p = ranksum(mean(PSTH_Cor_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t}),2)/binSize,mean(PSTH_Prem_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t}),2)/binSize);
                Sign.ITI_All.(trialTypo2{t})(i) = p;
%             end
            CorFR.All.(trialTypo2{t})(i) = mean(mean(PSTH_Cor_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t})/binSize,2)-mean(PSTH_Cor_All{1, VAR_ITIixOIMatched(i)}(:,1:10)/binSize,2));
            PremFR.All.(trialTypo2{t})(i) = mean(mean(PSTH_Prem_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t})/binSize,2)-mean(PSTH_Prem_All{1, VAR_ITIixOIMatched(i)}(:,1:10)/binSize,2));
            CorFR_AUC.All.(trialTypo2{t})(i,:) = sum(mean(PSTH_Cor_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t})/binSize,1)-mean(mean(PSTH_Cor_All{1, VAR_ITIixOIMatched(i)}(:,1:10)/binSize,1)));
            PremFR_AUC.All.(trialTypo2{t})(i,:) = sum(mean(PSTH_Prem_All{1, VAR_ITIixOIMatched(i)}(:,timeBins{t})/binSize,1)-mean(mean(PSTH_Prem_All{1, VAR_ITIixOIMatched(i)}(:,1:10)/binSize,1)));
    end
end
%Sign.ITI_All =  mafdr([Sign.ITI_All],'BHFDR',true);
%Sign.ITI_OI =  mafdr([Sign.ITI_OI],'BHFDR',true);


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
%Chi_RespTypesOIvsOther.pval = mafdr(Chi_RespTypesOIvsOther.pval,'BHFDR',true);


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
            %subplot(4,3,pos(n))
            hold on
            %plot([CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))); PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))],'k')
            %plot([CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)& Sign.ITI_OI.(trialTypo2{t}) < 0.05)); PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl) & Sign.ITI_OI.(trialTypo2{t}) < 0.05))],'r')
            %boxplot([CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))' PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))'],'Labels',{'Correct','Premature'})
            diffCorPrem.OI.(trialTypo2{t})(Cl) = median(((PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))'-CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))')./abs((CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))')))*100);
            %ylabel('Firing rate (Hz)')
            %diffCorPrem.OI2.(trialTypo2{t}){Cl} = ((PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))'-CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))')./(CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl)))'))*100;

            try
            CorVsPremSign.OI.(trialTypo2{t})(Cl) = signrank(CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))));
            %CorVsPremSign.OI.(trialTypo2{t})(Cl) = signrank(CorFR_AUC.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))),PremFR_AUC.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'==ClustUpMidDown(Cl))));

            %CorVsPremSign.OIVal.(trialTypo2{t}){Cl} = abs(CorFR.OI.(trialTypo2{t})(find(Sign.OI.(trialTypo2{t}) < 0.05))-PremFR.OI.(trialTypo2{t})(find(Sign.OI.(trialTypo2{t}) < 0.05)));
            catch
                CorVsPremSign.OI.(trialTypo2{t})(Cl) = nan;
            end
        end
        subplot(3,1,pos(n))
        bar(diffCorPrem.OI.(trialTypo2{t}))
        ylim([-25 25])
        n = n+1;
        %CorVsPremSign.OI.(trialTypo2{t}) = mafdr([CorVsPremSign.OI.(trialTypo2{t})],'BHFDR',true);
        [CorVsPremSign_Combined.OI.(trialTypo2{t}), ~, CorVsPremSign_Combined.OI.stats.(trialTypo2{t})] = signrank(abs(CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))),abs(PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2)))));
        CorVsPremSign_Combined2.OI.(trialTypo2{t}){2} = abs(CorFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2))));
        CorVsPremSign_Combined2.OI.(trialTypo2{t}){3} = abs(PremFR.OI.(trialTypo2{t})(find(ITI_OI_clust.(trialTypo2{t})'~=ClustUpMidDown(2))));
end
CorVsPremSign_Combined.OI_combined = [CorVsPremSign_Combined.OI.TH CorVsPremSign_Combined.OI.Resp_bef];
%CorVsPremSign_Combined.OI_combined = mafdr([CorVsPremSign_Combined.OI.TH CorVsPremSign_Combined.OI.Resp_bef],'BHFDR',true);
CorVsPremSign_Combined.OI.TH = CorVsPremSign_Combined.OI_combined(1);
CorVsPremSign_Combined.OI.Resp_bef = CorVsPremSign_Combined.OI_combined(2);
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

figure
n = 1;
pos = [1 2 3];
for t = [1 2 3]
        for Cl = 1:numel(ClustUpMidDown)
            %subplot(4,3,pos(n))
            %hold on
            %plot([CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))); PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))],'k')
            %plot([CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)& Sign.ITI_All.(trialTypo2{t}) < 0.05)); PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl) & Sign.ITI_All.(trialTypo2{t}) < 0.05))],'r')
            %boxplot([CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))' PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))'],'Labels',{'Correct','Premature'})
            %ylabel('Firing rate (Hz)')
            diffCorPrem.All.(trialTypo2{t})(Cl) = median(((PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))'-CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))')./abs((CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))')))*100);
            %diffCorPrem.All2.(trialTypo2{t}){Cl} = ((PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))'-CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))')./(CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))'))*100;

            %n = n+1;
            h1(i) = lillietest(CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))));
            h2(i) = lillietest(PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl)))); 
            if h1(i) == 0 & h2(i) == 0
                try
                [h CorVsPremSign.All.(trialTypo2{t})(Cl)] = ttest(CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))),PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))),'Vartype','unequal');
                %CorVsPremSign.AllVal.(trialTypo2{t}){Cl} = abs(CorFR.All.(trialTypo2{t})(find(Sign.ITI_All.(trialTypo2{t}) < 0.05))-PremFR.All.(trialTypo2{t})(find(Sign.ITI_All.(trialTypo2{t}) < 0.05)));

                catch
                    CorVsPremSign.All.(trialTypo2{t})(Cl) = nan;
                end
            else
                try
                CorVsPremSign.All.(trialTypo2{t})(Cl) = signrank(CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))),PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'==ClustUpMidDown(Cl))));
                %CorVsPremSign.AllVal.(trialTypo2{t}){Cl} = abs(CorFR.All.(trialTypo2{t})(find(Sign.ITI_All.(trialTypo2{t}) < 0.05))-PremFR.All.(trialTypo2{t})(find(Sign.ITI_All.(trialTypo2{t}) < 0.05)));

                catch
                    CorVsPremSign.All.(trialTypo2{t})(Cl) = nan;
            end
        end
        end
        subplot(3,1,pos(n))
        bar(diffCorPrem.All.(trialTypo2{t}))
        ylim([-25 25])
        n = n+1;
        %CorVsPremSign.All.(trialTypo2{t}) = mafdr([CorVsPremSign.All.(trialTypo2{t})],'BHFDR',true);
        CorVsPremSign_Combined.All.(trialTypo2{t}) = signrank(abs(CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'~=ClustUpMidDown(2)))),abs(PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'~=ClustUpMidDown(2)))));
        CorVsPremSign_Combined2.All.(trialTypo2{t}){2} = abs(CorFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'~=ClustUpMidDown(2))));
        CorVsPremSign_Combined2.All.(trialTypo2{t}){3} = abs(PremFR.All.(trialTypo2{t})(find(ITI_All_clust.(trialTypo2{t})(VAR_ITIixOIMatched)'~=ClustUpMidDown(2))));
end
CorVsPremSign_Combined.All_combined = mafdr([CorVsPremSign_Combined.All.TH CorVsPremSign_Combined.All.Resp_bef],'BHFDR',true);
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


%% 4. Get Z score normalized PSTHs for plotting colormap
%timeBins = {11:15, 11:20, 1:10, 12:20, 1:10};

        
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

figure
colormap('jet')
imagesc(B_OI(ib_OI(I),1:3))

figure
subplot(1,3,1)
colormap('jet')
%[B,I] = sort(mean(OI_PSTH.PSTH.TS.Cor(:,11:15),2),'descend');
imagesc(time_axis_TS,1:(size(OI_PSTH.PSTH.TS.Cor,1)),OI_PSTH.PSTH.TS.Cor(I,1:85))
axis([-2-(0.5*binSize) 15+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.TS.Cor,1)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Trial start')

% subplot(2,3,4)
% colormap('jet')
% imagesc(time_axis_TS,1:(size(OI_PSTH.PSTH.TS.Prem,1)),OI_PSTH.PSTH.TS.Prem(I,1:85))
% axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.TS.Prem,1)+0.5])
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
% imagesc(time_axis_TH,1:(size(OI_PSTH.PSTH.TH.Prem,1)),OI_PSTH.PSTH.TH.Prem(I,1:25))
% axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.TH.Prem,1)+0.5])
% ylabel('Neuron')
% caxis(axisRange)
% title('Cue orientation')

subplot(1,3,3)
colormap('jet')
%[B,I] = sort(mean(OI_PSTH.PSTH.Resp_bef.Cor(:,1:10),2),'descend');
imagesc(time_axis_Resp_bef,1:(size(OI_PSTH.PSTH.Resp_bef.Cor,1)),OI_PSTH.PSTH.Resp_bef.Cor(I,1:20))
axis([-2-(0.5*binSize) 2+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.Resp_bef.Cor,1)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Response')

% subplot(2,3,6)
% colormap('jet')
% imagesc(time_axis_Resp_bef,1:(size(OI_PSTH.PSTH.Resp_bef.Prem,1)),OI_PSTH.PSTH.Resp_bef.Prem(I,1:30))
% axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(OI_PSTH.PSTH.Resp_bef.Prem,1)+0.5])
% ylabel('Neuron')
% caxis(axisRange)
% title('Response')

% Plot Average PSTH per group
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


%[B,I] = sort(mean(All_PSTH.PSTH.TH.Cor(VAR_ITIixOIMatched,13:17),2),'descend');
I_oneIxs = find(ClustType_All(VAR_ITIixOIMatched)==1);
I_zeroIxs = find(ClustType_All(VAR_ITIixOIMatched)==0);
I_negIxs = find(ClustType_All(VAR_ITIixOIMatched)==-1);
[B,I_one] = sort(mean(All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched(I_oneIxs),1:5),2),'descend');
[B,I_zero] = sort(mean(All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched(I_zeroIxs),1:5),2),'descend');
[B,I_neg] = sort(mean(All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched(I_negIxs),1:5),2),'descend');
I = [I_oneIxs(I_one) I_zeroIxs(I_zero) I_negIxs(I_neg)];

figure
colormap('jet')
imagesc(B_All(ib_All(VAR_ITIixOIMatched(I)),1:3))

figure
subplot(1,3,1)
colormap('jet')
%[B,I] = sort(mean(All_PSTH.PSTH.TS.Cor(VAR_ITIixOIMatched,11:15),2),'descend');
imagesc(time_axis_TS,1:(size(All_PSTH.PSTH.TS.Cor(VAR_ITIixOIMatched),2)),All_PSTH.PSTH.TS.Cor(VAR_ITIixOIMatched(I),1:85))
axis([-2-(0.5*binSize) 15+(0.5*binSize) 0.5 size(All_PSTH.PSTH.TS.Cor(VAR_ITIixOIMatched),2)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Trial start')

% subplot(1,3,4)
% colormap('jet')
% imagesc(time_axis_TS,1:(size(All_PSTH.PSTH.TS.Prem(VAR_ITIixOIMatched),2)),All_PSTH.PSTH.TS.Prem(VAR_ITIixOIMatched(I),1:85))
% axis([-2-(0.5*binSize) 15+(0.5*binSize) 0.5 size(All_PSTH.PSTH.TS.Prem(VAR_ITIixOIMatched),2)+0.5])
% ylabel('Neuron')
% caxis(axisRange)
% title('Trial start')

subplot(1,3,2)
colormap('jet')
%[B,I] = sort(mean(All_PSTH.PSTH.TH.Cor(VAR_ITIixOIMatched,11:15),2),'descend');
imagesc(time_axis_TH,1:(size(All_PSTH.PSTH.TH.Cor(VAR_ITIixOIMatched),2)),All_PSTH.PSTH.TH.Cor(VAR_ITIixOIMatched(I),1:25))
axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(All_PSTH.PSTH.TH.Cor(VAR_ITIixOIMatched),2)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Cue orientation')

% subplot(2,3,5)
% colormap('jet')
% imagesc(time_axis_TH,1:(size(All_PSTH.PSTH.TH.Prem(VAR_ITIixOIMatched),2)),All_PSTH.PSTH.TH.Prem(VAR_ITIixOIMatched(I),1:25))
% axis([-2-(0.5*binSize) 3+(0.5*binSize) 0.5 size(All_PSTH.PSTH.TH.Prem(VAR_ITIixOIMatched),2)+0.5])
% ylabel('Neuron')
% caxis(axisRange)
% title('Cue orientation')

subplot(1,3,3)
colormap('jet')
%[B,I] = sort(mean(All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched,1:10),2),'descend');
imagesc(time_axis_Resp_bef,1:(size(All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched),2)),All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched(I),1:20))
axis([-2-(0.5*binSize) 2+(0.5*binSize) 0.5 size(All_PSTH.PSTH.Resp_bef.Cor(VAR_ITIixOIMatched),2)+0.5])
ylabel('Neuron')
caxis(axisRange)
title('Response')

% subplot(2,3,6)
% colormap('jet')
% imagesc(time_axis_Resp_bef,1:(size(All_PSTH.PSTH.Resp_bef.Prem(VAR_ITIixOIMatched),2)),All_PSTH.PSTH.Resp_bef.Prem(VAR_ITIixOIMatched(I),1:30))
% axis([-2-(0.5*binSize) 4+(0.5*binSize) 0.5 size(All_PSTH.PSTH.Resp_bef.Prem(VAR_ITIixOIMatched),2)+0.5])
% ylabel('Neuron')
% caxis(axisRange)
% title('Response')

%% Make plot that compares proportion of Cor vs Prem modualted units per sync point for OI and All

Prop.TSMod_OI = numel(find(Sign.ITI_OI.TS<0.05))/(numel(Sign.ITI_OI.TS));
Prop.THMod_OI  = numel(find(Sign.ITI_OI.TH<0.05))/(numel(Sign.ITI_OI.TH));
Prop.Resp_befMod_OI  = numel(find(Sign.ITI_OI.Resp_bef<0.05))/(numel(Sign.ITI_OI.Resp_bef));
Prop.TSMod_All = numel(find(Sign.ITI_All.TS<0.05))/(numel(Sign.ITI_All.TS));
Prop.THMod_All  = numel(find(Sign.ITI_All.TH<0.05))/(numel(Sign.ITI_All.TH));
Prop.Resp_befMod_All  = numel(find(Sign.ITI_All.Resp_bef<0.05))/(numel(Sign.ITI_All.Resp_bef));
figure
bar([Prop.TSMod_OI Prop.TSMod_All;Prop.THMod_OI Prop.THMod_All; Prop.Resp_befMod_OI Prop.Resp_befMod_All])
title('Proportion of Cor-Prem modulated neurons')
legend({'Identified','other'})

%% Get cue orientation times and cue orientation to response times
allFiles = dir;

n = 1;
VAR_ITIix = [];
VAR_SDix = [];
for file = 1:size(allFiles,1)
    if ~isempty(strfind(allFiles(file).name, 'txt'))
        MED_PC{n,1} = allFiles(file).name;
        n = n+1;
    end
end

for i = 1:size(neuronIndex.PyrIx_VAR_ITI_OI,2)
    recDate_ITI{i} = All{5,1}{neuronIndex.PyrIx_VAR_ITI_OI(i),2};
end
[Ix_ITI C_ITI IC_ITI] = unique(recDate_ITI);

n = 1;
for sessionITI = 1:size(C_ITI,1)
    [ITI.Correct{n,1}, ITI.Incorrect{n,1}, ITI.Omission{n,1}, PrematureITI(n,1), ITI.Perserverative{n,1}, nrTrialsStageITI(n,1), Perf_ITI{n,1}]  = CC_ephys_var_ITI_SD_TTL_Master(10,1,recDate_ITI{C_ITI(sessionITI)});
    n = n+1;
end

n = 2;
TH_Time_prem{1} = ITI_OI.ITI125.TrialDist{1, 1}.PremTreshold;
TH_Time_cor{1} = ITI_OI.ITI125.TrialDist{1, 1}.CorTreshold;
for neuron = 2:size(C_ITI,1)
    TH_Time_prem{n} = ITI_OI.ITI125.TrialDist{1, C_ITI(neuron)}.PremTreshold;
    TH_Time_cor{n} = ITI_OI.ITI125.TrialDist{1, C_ITI(neuron)}.CorTreshold;
    n = n+1;
end

% Get cue orientation times
Lat_TH_prem = [];
Lat_TH_cor = [];
n = 1;
for i = 1:size(TH_Time_prem,2)
    Lat_TH_prem = [Lat_TH_prem; (TH_Time_prem{i})];
    Lat_TH_cor = [Lat_TH_cor; (TH_Time_cor{i})];
    n = n+1;
end

% Get cue orientation to response time
Lat_THtoResp_prem = [];
for i = 1:size(C_ITI,1)
    Lat_THtoResp_prem = [Lat_THtoResp_prem (PrematureITI(i).LatencyITI125)];
end

Lat_THtoResp_cor = [];
for neuron = 1:size(C_ITI,1)
    Lat_THtoResp_cor = [Lat_THtoResp_cor (All{22, 1}{neuronIndex.PyrIx_VAR_ITI_OI(C_ITI(neuron)), 1}.corITI125Resp(ITI_OI.ITI125.TrialDist{1, C_ITI(neuron)}.CorTresholdIncluded)-All{22, 1}{neuronIndex.PyrIx_VAR_ITI_OI(C_ITI(neuron)), 1}.corITI125Cue(ITI_OI.ITI125.TrialDist{1, C_ITI(neuron)}.CorTresholdIncluded))+12.5];
end





% 
% 
%     
%     %% Plot Prem resp times and Cue orientation times
% allFiles = dir;
% 
% n = 1;
% VAR_ITIix = [];
% VAR_SDix = [];
% for file = 1:size(allFiles,1)
%     if ~isempty(strfind(allFiles(file).name, 'txt'))
%         MED_PC{n,1} = allFiles(file).name;
%         n = n+1;
%     end
% end
% 
% % for txtFile = 1:numel(MED_PC)
% %     for i = 1:size(All{5, 1},1)
% %         if All{5, 1}{i,2} == MED_PC{txtFile,1} & contains(All{5, 1}{i,1}, 'ITI')
% %             VAR_ITIix(end+1) = txtFile;
% %         elseif All{5, 1}{i,2} == MED_PC{txtFile,1} & contains(All{5, 1}{i,1}, 'SD')
% %             VAR_SDix(end+1) = txtFile;
% %         end
% %     end
% % end
% % 
% % VAR_ITIix = unique(VAR_ITIix);
% 
% 
% % Get performance per session
% % n = 1;
% % for sessionITI = VAR_ITIix
% %     [ITI.Correct{n,1}, ITI.Incorrect{n,1}, ITI.Omission{n,1}, PrematureITI(n,1), ITI.Perserverative{n,1}, nrTrialsStageITI(n,1), Perf_ITI{n,1}]  = CC_ephys_var_ITI_SD_TTL_Master(10,1,MED_PC{sessionITI,1});
% %     n = n+1;
% % end
% 
% n = 1;
% for sessionITI = 1:size(Ix_SD,2)
%     [ITI.Correct{n,1}, ITI.Incorrect{n,1}, ITI.Omission{n,1}, PrematureITI(n,1), ITI.Perserverative{n,1}, nrTrialsStageITI(n,1), Perf_ITI{n,1}]  = CC_ephys_var_ITI_SD_TTL_Master(10,1,Ix_SD{sessionITI});
%     n = n+1;
% end
% 
% Lat.ITI50 = [];
% Lat.ITI75 = [];
% Lat.ITI125 = [];
% 
% 
% %% TEMP
% % Lat.ITI125 = [];
% % n = 1;
% % for i = 1:size(PrematureITI,1)
% %     %Lat.ITI50 = (PrematureITI(i).LatencyITI50);
% %     %Lat.ITI75 = (PrematureITI(i).LatencyITI75);
% %     Lat.ITI125 = [Lat.ITI125 (PrematureITI(i).LatencyITI125)];
% %     n = n+1;
% % end
% % 
% % Lat.CorITI125 = [];
% % n = 1;
% % for i = 1:size(PrematureITI,1)
% %     %Lat.ITI50 = (PrematureITI(i).LatencyITI50);
% %     %Lat.ITI75 = (PrematureITI(i).LatencyITI75);
% %     Lat.CorITI125 = [Lat.ITI125 (CorrectITI(i).LatencyITI125)];
% %     n = n+1;
% % end
% %%
% 
% 
% n = 1;
% for i = 1:size(PrematureITI,1)
%     Lat.ITI50(n,:) = histcounts((PrematureITI(i).LatencyITI50),0:0.5:12.5);
%     Lat.ITI75(n,:) = histcounts((PrematureITI(i).LatencyITI75),0:0.5:12.5);
%     Lat.ITI125(n,:) = histcounts((PrematureITI(i).LatencyITI125),0:0.5:12.5)/numel(PrematureITI(i).LatencyITI125);
%     n = n+1;
% end
% for i = 1:size(Lat.ITI50,1)
%     if Lat.ITI50(i,11) > 0
%         Lat.ITI50(i,10) = Lat.ITI50(i,10)+Lat.ITI50(i,11);
%         Lat.ITI50(i,11) = 0;
%     elseif Lat.ITI75(i,16) > 0
%         Lat.ITI75(i,15) = Lat.ITI75(i,15)+Lat.ITI50(i,16);
%         Lat.ITI75(i,16) = 0;
%     end
% end
% 
% n = 2;
% TH_Time_prem{1} = ITI_OI.ITI125.TrialDist{1, 1}.PremTreshold;
% TH_Time_cor{1} = ITI_OI.ITI125.TrialDist{1, 1}.CorTreshold;
% for neuron = 2:size(ITI_OI.ITI125.TrialDist,2)
%     if numel(ITI_OI.ITI125.TrialDist{1, neuron}.PremTreshold)~=numel(TH_Time_prem{n-1})
%         TH_Time_prem{n} = ITI_OI.ITI125.TrialDist{1, neuron}.PremTreshold;
%         TH_Time_cor{n} = ITI_OI.ITI125.TrialDist{1, neuron}.CorTreshold;
%         n = n+1;
%     end
% end
% 
%% TEMP
% n = 2;
% TH_Time_prem{1} = ITI_OI.ITI125.TrialDist{1, 1}.PremTreshold;
% TH_Time_cor{1} = ITI_OI.ITI125.TrialDist{1, 1}.CorTreshold;
% for neuron = 2:size(C_ITI,1)
%     TH_Time_prem{n} = ITI_OI.ITI125.TrialDist{1, C_ITI(neuron)}.PremTreshold;
%     TH_Time_cor{n} = ITI_OI.ITI125.TrialDist{1, C_ITI(neuron)}.CorTreshold;
%     n = n+1;
% end
% 
% Lat_TH_prem = [];
% Lat_TH_cor = [];
% n = 1;
% for i = 1:size(TH_Time_prem,2)
%     Lat_TH_prem = [Lat_TH_prem; (TH_Time_prem{i})];
%     Lat_TH_cor = [Lat_TH_cor; (TH_Time_cor{i})];
%     n = n+1;
% end
% %%
% 
% n = 1;
% for i = 1:size(TH_Time_prem,2)
%     Lat_TH_prem(n,:) = histcounts((TH_Time_prem{i}),0:0.5:12.5)/numel(TH_Time_prem{i});
%     Lat_TH_cor(n,:) = histcounts((TH_Time_cor{i}),0:0.5:12.5)/numel(TH_Time_cor{i});
%     n = n+1;
% end
% 
% n = 1;
% for i = 1:size(TH_Time_prem,2)
%     Lat_TH_prem2(n) = median(TH_Time_prem{i});
%     Lat_TH_cor2(n) = median(TH_Time_cor{i});
%     n = n+1;
% end
% 
% x_ITI125 = 0.5:0.5:12.5;
% 
% figure
% plot(x_ITI125,[median(Lat.ITI125,1)'])
% hold on
% plot(x_ITI125,[median(Lat_TH_prem,1)'])
% plot(x_ITI125,[median(Lat_TH_cor,1)'])
% title('Premature response time Var. ITI')
% xlabel('Time (sec)')
% ylabel('Responses/session (mean)')
% legend({'ITI 5.0'})
% 
% figure
% errorbar(x_ITI125,[mean(Lat_TH_cor,1)'],std(Lat_TH_cor,1))
% title('Premature response time Var. ITI')
% xlabel('Time (sec)')
% ylabel('Responses/session (mean)')
% legend({'ITI 5.0'})
% 
% figure
% hold on
% errorbar(x_ITI125,[mean(Lat_TH_cor,1)'],std(Lat_TH_cor,1))
% errorbar(x_ITI125,[mean(Lat_TH_prem,1)'],std(Lat_TH_prem,1))
% title('Premature response time Var. ITI')
% xlabel('Time (sec)')
% ylabel('Responses/session (mean)')
% %legend({'ITI 5.0'})
% 
% figure
% hold on
% boxplot([Lat_TH_cor2' Lat_TH_prem2'],'Labels',{'Correct','Premature'})
% plot([Lat_TH_cor2; Lat_TH_prem2])
% ylabel('Start trial - Cue orientation (sec)')
% %ylim([0 20])
% %scatter(ones(size(Lat_TH_cor2,1),1).*(1+(rand(size(Lat_TH_cor2,1),1)-0.5)/5),Lat_TH_cor2,20,'k','filled')
% %scatter(ones(size(Lat_TH_prem2,1),1).*(2+(rand(size(Lat_TH_prem2,1),1)-0.5)/5),Lat_TH_prem2,20,'k','filled')
% 
% % Paired stats on time to cue orientation
% [p_Prem_timetoCueOri,h_Prem_timetoCueOri,stats_Prem_timetoCueOri] = signrank(Lat_TH_cor2, Lat_TH_prem2);


%% Junk

% for t = [1 2 5]
%     tot_All_OI = 
%     tot_All_All
%     for Cl = 1:numel(ClustUpMidDown)
%         x = [zeros(tot_All_OI,1); ones(tot_All_All,1)]';
%         y = [zeros(numel(find([ITI_OI_clust] == ClustUpMidDown(i))),1); ones(tot_All_OI-numel(find([ITI_OI_clust] == ClustUpMidDown(i))),1); zeros(numel(find([ITI_All_clust] == ClustUpMidDown(i))),1); ones(tot_All_All-numel(find([ITI_All_clust] == ClustUpMidDown(i))),1)]';
%         [tbl_All{i},chi2stat_All{i},Chi_All.pval(i)] = crosstab(x,y);
%         if numel(find([ITI_OI_clust] == ClustUpMidDown(i))) < 5 || numel(find([ITI_All_clust] == ClustUpMidDown(i))) < 5
%             tbl_All{i} = table([numel(find([ITI_OI_clust] == ClustUpMidDown(i)));tot_All_OI-numel(find([ITI_OI_clust] == ClustUpMidDown(i)))],[numel(find([ITI_All_clust] == ClustUpMidDown(i)));tot_All_All-numel(find([ITI_All_clust] == ClustUpMidDown(i)))],'VariableNames',{'Flu','NoFlu'},'RowNames',{'NoShot','Shot'});
%             [~,Chi_All.pval(i),chi2stat_All{i}] = fishertest(tbl_All{i});
%         end
%     end
% Chi_All.pval = mafdr(Chi_All.pval,'BHFDR',true);
% 
% 
% 
% %% premature modulation index
% for i = 1:size(CorFR.OI.Resp_bef,2)
%     AMI_OI.Resp_bef(i) = (CorFR.OI.Resp_bef(i)-PremFR.OI.Resp_bef(i))/(CorFR.OI.Resp_bef(i)+PremFR.OI.Resp_bef(i));
% end
% for i = 1:size(CorFR.All.Resp_bef,2)
%     AMI_All.Resp_bef(i) = (CorFR.All.Resp_bef(i)-PremFR.All.Resp_bef(i))/(CorFR.All.Resp_bef(i)+PremFR.All.Resp_bef(i));
% end
% 
% for n = 1:size(ITI_OI.All.PSTH,2)
%     All_PSTH.TS.Cor(n,:) = (ITI_OI.(ITIperiod).PSTH{1, n}.Cor.trial_start-mean(ITI_OI.All.PSTH{1, n}.Cor.trial_start(1:10)))/max(abs(ITI_OI.All.PSTH{1, n}.Cor.trial_start-mean(ITI_OI.All.PSTH{1, n}.Cor.trial_start(1:10))));%,'gaussian',3)/SD_OI.All.ID{1, n}.Cor.TrialStart.Peak_real;
%     All_PSTH.TH.Cor(n,:) = (ITI_OI.(ITIperiod).PSTH{1, n}.Cor.wait_start-mean(ITI_OI.All.PSTH{1, n}.Cor.trial_start(1:10)))/max(abs(ITI_OI.All.PSTH{1, n}.Cor.wait_start-mean(ITI_OI.All.PSTH{1, n}.Cor.trial_start(1:10))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Treshold.Peak_real;
%     All_PSTH.Resp.Cor(n,:) = (ITI_OI.(ITIperiod).PSTH{1, n}.Cor.resp-mean(ITI_OI.All.PSTH{1, n}.Cor.trial_start(1:10)))/max(abs(ITI_OI.All.PSTH{1, n}.Cor.resp-mean(ITI_OI.All.PSTH{1, n}.Cor.trial_start(1:10))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Response.Peak_real;
%     All_PSTH.TS.Prem(n,:) = (ITI_OI.(ITIperiod).PSTH{1, n}.Prem.trial_start-mean(ITI_OI.All.PSTH{1, n}.Cor.trial_start(1:10)))/max(abs(ITI_OI.All.PSTH{1, n}.Cor.trial_start-mean(ITI_OI.All.PSTH{1, n}.Cor.trial_start(1:10))));%,'gaussian',3)/SD_OI.All.ID{1, n}.Cor.TrialStart.Peak_real;
%     All_PSTH.TH.Prem(n,:) = (ITI_OI.(ITIperiod).PSTH{1, n}.Prem.wait_start-mean(ITI_OI.All.PSTH{1, n}.Cor.trial_start(1:10)))/max(abs(ITI_OI.All.PSTH{1, n}.Cor.wait_start-mean(ITI_OI.All.PSTH{1, n}.Cor.trial_start(1:10))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Treshold.Peak_real;
%     All_PSTH.Resp.Prem(n,:) = (ITI_OI.(ITIperiod).PSTH{1, n}.Prem.resp-mean(ITI_OI.All.PSTH{1, n}.Cor.trial_start(1:10)))/max(abs(ITI_OI.All.PSTH{1, n}.Cor.resp-mean(ITI_OI.All.PSTH{1, n}.Cor.trial_start(1:10))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Response.Peak_real;
% end
% 
% for Cl = 1:numel(ClustUpMidDown)
% subplot(2,3,Cl+3)
% time_x = -2:binSize:12.8;
% plot(time_x,[mean(All_PSTH.TS.Cor(find(Sign.ITI_OI(1,:)<0.05&ClustType_OI(1,:)==ClustUpMidDown(Cl)),1:75),1)],'k')
% hold on
% plot(time_x,[mean(All_PSTH.TS.Prem(find(Sign.ITI_OI(1,:)<0.05&ClustType_OI(1,:)==ClustUpMidDown(Cl)),1:75),1)],'b')
% legend({'Cor','Prem'})
% xlabel('Time cue ori')
% ylabel('FR (cor max normalized)')
% end