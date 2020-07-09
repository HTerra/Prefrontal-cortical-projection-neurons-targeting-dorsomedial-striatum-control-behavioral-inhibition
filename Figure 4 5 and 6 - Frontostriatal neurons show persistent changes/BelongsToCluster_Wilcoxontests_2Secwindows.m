function [ClustType, All, Clust_ord] = BelongsToCluster_Wilcoxontests_2Secwindows(All, ITI_All, ITI_OI, neuronIndex, parameters_OI, SD_All, SD_OI)

% Function calculates if neurons are acivated, silenced or unchanged in time windows, after trial start, after cue orientation and before the cue.
% Significance is calculated with a Wilcoxon signed-rank test compared to a 2 sec baseline period.


% Does not exclude neurons with buddy. Does exclude sessions without OI
% neurons.
ITIperiod = 'ITI125';

% Get correct PETHs from all pyramidal neurons, OI and non OI
PSTH_All = [SD_All.All.PSTH SD_OI.All.PSTH ITI_All.(ITIperiod).PSTH ITI_OI.(ITIperiod).PSTH];
PSTHBinVal_All = [SD_All.All.psthBinsValue SD_OI.All.psthBinsValue ITI_All.(ITIperiod).psthBinsValue ITI_OI.(ITIperiod).psthBinsValue];
        
% Get baseline subtracted PETHs
for n = 1:size(PSTH_All,2)
    try
        All_PSTH.PSTH.TS.Cor(n,:) = (PSTH_All{1, n}.Cor.trial_start);%-mean(PSTH_All{1, n}.Cor.trial_start(1:10)));%/mean(std(PSTHBinVal_All{1, n}.Cor.trial_start(:,1:10),1),2);%max(abs(SD_OI.All.PSTH{1, n}.Cor.trial_start-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:5))));%,'gaussian',3)/SD_OI.All.ID{1, n}.Cor.TrialStart.Peak_real;
    catch
         All_PSTH.PSTH.TS.Cor(n,:) = (PSTH_All{1, n}.Cor.trial_start(1:45));
    end
    try
        All_PSTH.PSTH.TH.Cor(n,:) = (PSTH_All{1, n}.Cor.wait_start);%-mean(PSTH_All{1, n}.Cor.trial_start(1:10)));%/mean(std(PSTHBinVal_All{1, n}.Cor.trial_start(:,1:10),1),2);%max(abs(PSTH_SD_All.PSTH{1, n}.Cor.wait_start-mean(PSTH_SD_All.PSTH{1, n}.Cor.trial_start(1:5))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Treshold.Peak_real;
    catch
        All_PSTH.PSTH.TH.Cor(n,:) = (PSTH_All{1, n}.Cor.wait_start(1:35));%
    end
    All_PSTH.PSTH.Resp.Cor(n,:) = (PSTH_All{1, n}.Cor.resp);%-mean(PSTH_All{1, n}.Cor.trial_start(1:10)));%/mean(std(PSTHBinVal_All{1, n}.Cor.trial_start(:,1:10),1),2);%max(abs(PSTH_SD_All.PSTH{1, n}.Cor.resp-mean(PSTH_SD_All.PSTH{1, n}.Cor.trial_start(1:5))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Response.Peak_real;
    try
        All_PSTH.PSTH.Cue.Cor(n,:) = (PSTH_All{1, n}.Cor.cue);%-mean(PSTH_All{1, n}.Cor.trial_start(1:10)));%/mean(std(PSTHBinVal_All{1, n}.Cor.trial_start(:,1:10),1),2);%max(abs(PSTH_SD_All.PSTH{1, n}.Cor.resp-mean(PSTH_SD_All.PSTH{1, n}.Cor.trial_start(1:5))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Response.Peak_real;
    catch
        All_PSTH.PSTH.Cue.Cor(n,:) = (PSTH_All{1, n}.Cor.trial_start(26:45));%-mean(PSTH_All{1, n}.Cor.trial_start(1:10)));%/mean(std(PSTHBinVal_All{1, n}.Cor.trial_start(:,1:10),1),2);%max(abs(PSTH_SD_All.PSTH{1, n}.Cor.resp-mean(PSTH_SD_All.PSTH{1, n}.Cor.trial_start(1:5))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Response.Peak_real;
    end
end

for n = 1:size(PSTH_All,2)
    try
        All_PSTH.psthBinsValue.TS.Cor{n,:} = (PSTHBinVal_All{1, n}.Cor.trial_start);%/mean(std(PSTHBinVal_All{1, n}.Cor.trial_start(:,1:10),1),2);%max(abs(SD_OI.All.PSTH{1, n}.Cor.trial_start-mean(SD_OI.All.PSTH{1, n}.Cor.trial_start(1:5))));%,'gaussian',3)/SD_OI.All.ID{1, n}.Cor.TrialStart.Peak_real;
    catch
        All_PSTH.psthBinsValue.TS.Cor{n,:} = (PSTHBinVal_All{1, n}.Cor.trial_start(1:45));%
    end
    try
        All_PSTH.psthBinsValue.TH.Cor{n,:} = (PSTHBinVal_All{1, n}.Cor.wait_start);%/mean(std(PSTHBinVal_All{1, n}.Cor.trial_start(:,1:10),1),2);%max(abs(PSTH_SD_All.PSTH{1, n}.Cor.wait_start-mean(PSTH_SD_All.PSTH{1, n}.Cor.trial_start(1:5))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Treshold.Peak_real;
    catch
        All_PSTH.psthBinsValue.TH.Cor{n,:} = (PSTHBinVal_All{1, n}.Cor.wait_start(1:35));%
    end
        All_PSTH.psthBinsValue.Resp.Cor{n,:} = (PSTHBinVal_All{1, n}.Cor.resp);%/mean(std(PSTHBinVal_All{1, n}.Cor.trial_start(:,1:10),1),2);%max(abs(PSTH_SD_All.PSTH{1, n}.Cor.resp-mean(PSTH_SD_All.PSTH{1, n}.Cor.trial_start(1:5))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Response.Peak_real;
    try
        All_PSTH.psthBinsValue.Cue.Cor{n,:} = (PSTHBinVal_All{1, n}.Cor.cue);%/mean(std(PSTHBinVal_All{1, n}.Cor.trial_start(:,1:10),1),2);%max(abs(PSTH_SD_All.PSTH{1, n}.Cor.resp-mean(PSTH_SD_All.PSTH{1, n}.Cor.trial_start(1:5))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Response.Peak_real;
    catch
        All_PSTH.psthBinsValue.Cue.Cor{n,:} = (PSTHBinVal_All{1, n}.Cor.trial_start(:,26:45));%/mean(std(PSTHBinVal_All{1, n}.Cor.trial_start(:,1:10),1),2);%max(abs(PSTH_SD_All.PSTH{1, n}.Cor.resp-mean(PSTH_SD_All.PSTH{1, n}.Cor.trial_start(1:5))));%,'gaussian',3)/PSTH_SD_All.ID{1, n}.Cor.Response.Peak_real;
    end
end

Type.TS = [];
Type.Cue = [];
Type.Resp = [];
trialTypo = {'TS','TH','Cue','Resp','Resp'};
trialTypo2 = {'TS','TH','Cue','Resp','Resp_bef'};
timeBins = {11:15, 13:22, 1:10, 12:16, 1:10};
for t = 1:size(trialTypo,2)
    timeRange = timeBins{t};
    for i = 1:size(All_PSTH.psthBinsValue.(trialTypo{t}).Cor,1)
            h1(i) = lillietest(mean(All_PSTH.psthBinsValue.TS.Cor{i, 1}(:,1:10)/0.2,2));
            h2(i) = lillietest(mean(All_PSTH.psthBinsValue.(trialTypo{t}).Cor{i, 1}(:,timeRange)/0.2,2)); 

            if h1(i) == 0 & h2(i) == 0
                [~, p] = ttest(mean(All_PSTH.psthBinsValue.TS.Cor{i, 1}(:,1:10)/0.2,2),mean(All_PSTH.psthBinsValue.(trialTypo{t}).Cor{i, 1}(:,timeRange)/0.2,2));
                Type.(trialTypo2{t})(i,1) = p;
                Type.(trialTypo2{t})(i,2) = mean(mean(All_PSTH.psthBinsValue.(trialTypo{t}).Cor{i, 1}(:,timeRange)/0.2,2))-mean(mean(All_PSTH.psthBinsValue.TS.Cor{i, 1}(:,1:10)/0.2,2));
            else
                p = signrank(mean(All_PSTH.psthBinsValue.TS.Cor{i, 1}(:,1:10)/0.2,2),mean(All_PSTH.psthBinsValue.(trialTypo{t}).Cor{i, 1}(:,timeRange)/0.2,2));
                Type.(trialTypo2{t})(i,1) = p;
                Type.(trialTypo2{t})(i,2) = median(mean(All_PSTH.psthBinsValue.(trialTypo{t}).Cor{i, 1}(:,timeRange)/0.2,2))-median(mean(All_PSTH.psthBinsValue.TS.Cor{i, 1}(:,1:10)/0.2,2));

            end
    end
    for i = 1:size(All_PSTH.psthBinsValue.(trialTypo{t}).Cor,1)
        if Type.(trialTypo2{t})(i,1) < 0.05 & Type.(trialTypo2{t})(i,2) > 0
            Type.(trialTypo2{t})(i,1) = 1;
        elseif Type.(trialTypo2{t})(i,1) < 0.05 & Type.(trialTypo2{t})(i,2) < 0
            Type.(trialTypo2{t})(i,1) = -1;
        else
            Type.(trialTypo2{t})(i,1) = 0;
        end
    end
end

figure
bar([numel(find(Type.TS(:,1)==1)) numel(find(Type.TS(:,1)==0)) numel(find(Type.TS(:,1)==-1)); numel(find(Type.TH(:,1)==1)) numel(find(Type.TH(:,1)==0)) numel(find(Type.TH(:,1)==-1)); numel(find(Type.Cue(:,1)==1)) numel(find(Type.Cue(:,1)==0)) numel(find(Type.Cue(:,1)==-1)); numel(find(Type.Resp(:,1)==1)) numel(find(Type.Resp(:,1)==0)) numel(find(Type.Resp(:,1)==-1))],'stacked');

for t = 1:size(trialTypo,2)
    % Get which neuron belongs to which cluster
    SD_All_clust.(trialTypo2{t}) = Type.(trialTypo2{t})(1:size(SD_All.All.PSTH,2),1);
    SD_OI_clust.(trialTypo2{t}) = Type.(trialTypo2{t})(size(SD_All.All.PSTH,2)+1:size(SD_All.All.PSTH,2)+size(SD_OI.All.PSTH,2),1);
    ITI_All_clust.(trialTypo2{t}) = Type.(trialTypo2{t})(size(SD_All.All.PSTH,2)+size(SD_OI.All.PSTH,2)+1:size(SD_All.All.PSTH,2)+size(SD_OI.All.PSTH,2)+size(neuronIndex.PyrIxVAR_ITI  ,2),1);
    ITI_OI_clust.(trialTypo2{t}) = Type.(trialTypo2{t})(size(SD_All.All.PSTH,2)+size(SD_OI.All.PSTH,2)+size(neuronIndex.PyrIxVAR_ITI  ,2)+1:size(SD_All.All.PSTH,2)+size(SD_OI.All.PSTH,2)+size(neuronIndex.PyrIxVAR_ITI  ,2)+size(neuronIndex.PyrIx_VAR_ITI_OI,2),1);
    
    SD_All_clust_deltaFR.(trialTypo2{t}) = Type.(trialTypo2{t})(1:size(SD_All.All.PSTH,2),2);
    SD_OI_clust_deltaFR.(trialTypo2{t}) = Type.(trialTypo2{t})(size(SD_All.All.PSTH,2)+1:size(SD_All.All.PSTH,2)+size(SD_OI.All.PSTH,2),2);
    ITI_All_clust_deltaFR.(trialTypo2{t}) = Type.(trialTypo2{t})(size(SD_All.All.PSTH,2)+size(SD_OI.All.PSTH,2)+1:size(SD_All.All.PSTH,2)+size(SD_OI.All.PSTH,2)+size(neuronIndex.PyrIxVAR_ITI  ,2),2);
    ITI_OI_clust_deltaFR.(trialTypo2{t}) = Type.(trialTypo2{t})(size(SD_All.All.PSTH,2)+size(SD_OI.All.PSTH,2)+size(neuronIndex.PyrIxVAR_ITI  ,2)+1:size(SD_All.All.PSTH,2)+size(SD_OI.All.PSTH,2)+size(neuronIndex.PyrIxVAR_ITI  ,2)+size(neuronIndex.PyrIx_VAR_ITI_OI,2),2);


    % Make struct for function output
    ClustType.(trialTypo2{t}).SD_All_clust = SD_All_clust;
    ClustType.(trialTypo2{t}).SD_OI_clust = SD_OI_clust;
    ClustType.(trialTypo2{t}).ITI_All_clust = ITI_All_clust;
    ClustType.(trialTypo2{t}).ITI_OI_clust = ITI_OI_clust;

    % Fill in cluster type in All struct
    %All{24,1} = nan;
    All{28,1}(neuronIndex.PyrIx_VAR_ITI_OI,t+1) = ITI_OI_clust.(trialTypo2{t});
    All{28,1}(neuronIndex.PyrIxVAR_ITI,t+1) = ITI_All_clust.(trialTypo2{t});
    All{28,1}(neuronIndex.PyrIx_VAR_SD_OI,t+1) = SD_OI_clust.(trialTypo2{t});
    All{28,1}(neuronIndex.PyrIxVAR_SD,t+1) = SD_All_clust.(trialTypo2{t});
    All{28,2} ='twoSecWindow';
end

VAR_SDixOIMatched = [];
for neuron = 1:numel(neuronIndex.PyrIxVAR_SD)
    for i = neuronIndex.PyrIx_VAR_SD_OI
        if strcmp(All{9, 1}{neuronIndex.PyrIxVAR_SD(neuron),1},All{9, 1}{i,1})
            VAR_SDixOIMatched = [VAR_SDixOIMatched neuron];
            break
        else
        end
    end
end
VAR_ITIixOIMatched = [];
for neuron = 1:numel(neuronIndex.PyrIxVAR_ITI)
    for i = neuronIndex.PyrIx_VAR_ITI_OI
        if strcmp(All{9, 1}{neuronIndex.PyrIxVAR_ITI(neuron),1},All{9, 1}{i,1})
            VAR_ITIixOIMatched = [VAR_ITIixOIMatched neuron];
            break
        else
        end
    end
end

for t = 1:size(trialTypo,2)
    Clust_ord = [1 0 -1];
    tot_SD_OI = size(SD_OI_clust.(trialTypo{t}),1);
    tot_SD_All = size(SD_All_clust.(trialTypo{t})(VAR_SDixOIMatched),1);
    for i = 1:3
        x = [zeros(tot_SD_OI,1); ones(tot_SD_All,1)]';
        y = [zeros(numel(find(SD_OI_clust.(trialTypo{t}) == Clust_ord(i))),1); ones(tot_SD_OI-numel(find(SD_OI_clust.(trialTypo{t}) == Clust_ord(i))),1); zeros(numel(find(SD_All_clust.(trialTypo{t})(VAR_SDixOIMatched) == Clust_ord(i))),1); ones(tot_SD_All-numel(find(SD_All_clust.(trialTypo{t})(VAR_SDixOIMatched) == Clust_ord(i))),1)]';
        [tbl_SD.(trialTypo{t}){i},chi2stat_SD.(trialTypo{t}){i},Chi_SD.pval.(trialTypo{t})(i)] = crosstab(x,y);
        if numel(find(SD_OI_clust.(trialTypo{t}) == Clust_ord(i))) < 5 || numel(find(SD_All_clust.(trialTypo{t})(VAR_SDixOIMatched) == Clust_ord(i))) < 5
            tbl_SD.(trialTypo{t}){i} = table([numel(find(SD_OI_clust.(trialTypo{t}) == Clust_ord(i)));tot_SD_OI-numel(find(SD_OI_clust.(trialTypo{t}) == Clust_ord(i)))],[numel(find(SD_All_clust.(trialTypo{t})(VAR_SDixOIMatched) == Clust_ord(i)));tot_SD_All-numel(find(SD_All_clust.(trialTypo{t})(VAR_SDixOIMatched) == Clust_ord(i)))],'VariableNames',{'Flu','NoFlu'},'RowNames',{'NoShot','Shot'});
            [~,Chi_SD.pval.(trialTypo{t})(i),chi2stat_SD.(trialTypo{t}){i}] = fishertest(tbl_SD.(trialTypo{t}){i});
        end
    end
end
for t = 1:size(trialTypo,2)
    Clust_ord = [1 0 -1];
    tot_ITI_OI = size(ITI_OI_clust.(trialTypo{t}),1);
    tot_ITI_All = size(ITI_All_clust.(trialTypo{t})(VAR_ITIixOIMatched),1);
    for i = 1:3
        x = [zeros(tot_ITI_OI,1); ones(tot_ITI_All,1)]';
        y = [zeros(numel(find(ITI_OI_clust.(trialTypo{t}) == Clust_ord(i))),1); ones(tot_ITI_OI-numel(find(ITI_OI_clust.(trialTypo{t}) == Clust_ord(i))),1); zeros(numel(find(ITI_All_clust.(trialTypo{t})(VAR_ITIixOIMatched) == Clust_ord(i))),1); ones(tot_ITI_All-numel(find(ITI_All_clust.(trialTypo{t})(VAR_ITIixOIMatched) == Clust_ord(i))),1)]';
        [tbl_ITI.(trialTypo{t}){i},chi2stat_ITI.(trialTypo{t}){i},Chi_ITI.pval.(trialTypo{t})(i)] = crosstab(x,y);
        if numel(find(ITI_OI_clust.(trialTypo{t}) == Clust_ord(i))) < 5 || numel(find(ITI_All_clust.(trialTypo{t})(VAR_ITIixOIMatched) == Clust_ord(i))) < 5
            tbl_ITI.(trialTypo{t}){i} = table([numel(find(ITI_OI_clust.(trialTypo{t}) == Clust_ord(i)));tot_ITI_OI-numel(find(ITI_OI_clust.(trialTypo{t}) == Clust_ord(i)))],[numel(find(ITI_All_clust.(trialTypo{t})(VAR_ITIixOIMatched) == Clust_ord(i)));tot_ITI_All-numel(find(ITI_All_clust.(trialTypo{t})(VAR_ITIixOIMatched) == Clust_ord(i)))],'VariableNames',{'Flu','NoFlu'},'RowNames',{'NoShot','Shot'});
            [~,Chi_ITI.pval.(trialTypo{t})(i),chi2stat_ITI.(trialTypo{t}){i}] = fishertest(tbl_ITI.(trialTypo{t}){i});
        end
    end
end
