function [SD, ITI, SPSS, nrTrialsStageITI, nrTrialsStageSD, PrematureITI, PrematureSD, MED_PC] = Behavior_Analysis(All,neuronIndex)

% Find all .txt files (assumed to be MED_PC file) and sort them in VAR_SD
% and VAR_ITI MED_PC files for further analysis
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

for txtFile = 1:numel(MED_PC)
    for i = 1:size(All{5, 1},1)
        if All{5, 1}{i,2} == MED_PC{txtFile,1} & contains(All{5, 1}{i,1}, 'ITI')
            VAR_ITIix(end+1) = txtFile;
        elseif All{5, 1}{i,2} == MED_PC{txtFile,1} & contains(All{5, 1}{i,1}, 'SD')
            VAR_SDix(end+1) = txtFile;
        end
    end
end

VAR_ITIix = unique(VAR_ITIix);
VAR_SDix = unique(VAR_SDix);

% Get performance per session
n = 1;
for sessionITI = VAR_ITIix
    [ITI.Correct{n,1}, ITI.Incorrect{n,1}, ITI.Omission{n,1}, PrematureITI(n,1), ITI.Perserverative{n,1}, nrTrialsStageITI(n,1), Perf_ITI{n,1}]  = CC_ephys_var_ITI_SD_TTL_Master(10,1,MED_PC{sessionITI,1});
    n = n+1;
end

n = 1;
for sessionSD = VAR_SDix
    [SD.Correct{n,1}, SD.Incorrect{n,1}, SD.Omission{n,1}, PrematureSD(n,1), SD.Perserverative{n,1}, nrTrialsStageSD(n,1), Perf_SD{n,1}]  = CC_ephys_var_ITI_SD_TTL_Master(10,1,MED_PC{sessionSD,1});
    n = n+1;
end


%% Plot behavioral performance for VAR_SD
for i = 1:size(Perf_SD,1)
GemAc(i,1) = Perf_SD{i,1}.accuracyPerSD(1,1);
GemAc(i,2) = Perf_SD{i,1}.accuracyPerSD(1,2);
GemAc(i,3) = Perf_SD{i,1}.accuracyPerSD(1,3);

GemOm(i,1) = Perf_SD{i,1}.omissionPerSD(1,1);
GemOm(i,2) = Perf_SD{i,1}.omissionPerSD(1,2);
GemOm(i,3) = Perf_SD{i,1}.omissionPerSD(1,3);

GemPrem(i,1) = Perf_SD{i,1}.prematurePerSD(1,1);
GemPrem(i,2) = Perf_SD{i,1}.prematurePerSD(1,2);
GemPrem(i,3) = Perf_SD{i,1}.prematurePerSD(1,3);

GemCor(i,1) = numel(SD.Correct{i, 1}.TrialPosSD2);
GemCor(i,2) = numel(SD.Correct{i, 1}.TrialPosSD5);
GemCor(i,3) = numel(SD.Correct{i, 1}.TrialPosSD10);

GemInc(i,1) = numel(SD.Incorrect{i, 1}.TrialPosSD2);
GemInc(i,2) = numel(SD.Incorrect{i, 1}.TrialPosSD5);
GemInc(i,3) = numel(SD.Incorrect{i, 1}.TrialPosSD10);
end

% Find first position of unique pattern in variable ai
[~, ai, ~] = unique(GemOm(:,1),'stable');

for i = 1:size(nrTrialsStageSD,1)
    nrTrSD(i) = nrTrialsStageSD(i,1).TotalNr;
end


% Make boxplots and do kruskall wallis statistics with multiple comparisons
figure
subplot(2,4,1)
hold on
%notBoxPlot([GemAc(ai,1) GemAc(ai,2) GemAc(ai,3)])
plot([GemAc(ai,3)'; GemAc(ai,2)'; GemAc(ai,1)'],'k')
boxplot([GemAc(ai,3) GemAc(ai,2) GemAc(ai,1)],'Labels',{'1.0','0.5','0.20'})
%scatter(ones(size(GemAc(ai,1))).*(1+(rand(size(GemAc(ai,1)))-0.5)/5),GemAc(ai,1),20,'r','filled')
%scatter(ones(size(GemAc(ai,2))).*(2+(rand(size(GemAc(ai,2)))-0.5)/5),GemAc(ai,2),20,'r','filled')
%scatter(ones(size(GemAc(ai,3))).*(3+(rand(size(GemAc(ai,3)))-0.5)/5),GemAc(ai,3),20,'r','filled')
ylabel('Accuracy (%)')
xlabel('SD (sec)');
ylim([50 100])

% subplot(2,4,1)
% hold on
% %notBoxPlot([GemAc(ai,1) GemAc(ai,2) GemAc(ai,3)])
% plot([(GemCor(ai,3)./nrTrSD'*100)'; (GemCor(ai,2)./nrTrSD'*100)'; (GemCor(ai,1)./nrTrSD'*100)'],'k')
% boxplot([(GemCor(ai,3)./nrTrSD'*100) (GemCor(ai,2)./nrTrSD'*100) (GemCor(ai,1)./nrTrSD'*100)],'Labels',{'1.0','0.5','0.20'})
% %scatter(ones(size(GemAc(ai,1))).*(1+(rand(size(GemAc(ai,1)))-0.5)/5),GemAc(ai,1),20,'r','filled')
% %scatter(ones(size(GemAc(ai,2))).*(2+(rand(size(GemAc(ai,2)))-0.5)/5),GemAc(ai,2),20,'r','filled')
% %scatter(ones(size(GemAc(ai,3))).*(3+(rand(size(GemAc(ai,3)))-0.5)/5),GemAc(ai,3),20,'r','filled')
% ylabel('Correct (%)')
% xlabel('SD (sec)');
% ylim([50 100])

% subplot(2,4,2)
% hold on
% %notBoxPlot([GemAc(ai,1) GemAc(ai,2) GemAc(ai,3)])
% plot([(GemInc(ai,3)./nrTrSD'*100)'; (GemInc(ai,2)./nrTrSD'*100)'; (GemInc(ai,1)./nrTrSD'*100)'],'k')
% boxplot([(GemInc(ai,3)./nrTrSD'*100) (GemInc(ai,2)./nrTrSD'*100) (GemInc(ai,1)./nrTrSD'*100)],'Labels',{'1.0','0.5','0.20'})
% %scatter(ones(size(GemAc(ai,1))).*(1+(rand(size(GemAc(ai,1)))-0.5)/5),GemAc(ai,1),20,'r','filled')
% %scatter(ones(size(GemAc(ai,2))).*(2+(rand(size(GemAc(ai,2)))-0.5)/5),GemAc(ai,2),20,'r','filled')
% %scatter(ones(size(GemAc(ai,3))).*(3+(rand(size(GemAc(ai,3)))-0.5)/5),GemAc(ai,3),20,'r','filled')
% ylabel('Incorrect (%)')
% xlabel('SD (sec)');
% ylim([0 100])

SPSS.SD.Acc = [GemAc(ai,1), GemAc(ai,2), GemAc(ai,3)];
%GemAcStatGroup = [ones(size(GemAc(ai,1),1),1);  2*ones(size(GemAc(ai,2),1),1); 3*ones(size(GemAc(ai,3),1),1)];
%[AcStat{1},AcStat{2},AcStat{3}] = kruskalwallis(GemAcStat,GemAcStatGroup,'off');
%AcStatSD{4} = multcompare(AcStat{3},'Display','off');

subplot(2,4,2)
hold on
%notBoxPlot([GemOm(ai,1) GemOm(ai,2) GemOm(ai,3)])
plot([GemOm(ai,3)'; GemOm(ai,2)'; GemOm(ai,1)'],'k')
boxplot([GemOm(ai,3) GemOm(ai,2) GemOm(ai,1)],'Labels',{'1.0','0.5','0.20'})
% scatter(ones(size(GemOm(ai,1))).*(1+(rand(size(GemOm(ai,1)))-0.5)/5),GemOm(ai,1),20,'r','filled')
% scatter(ones(size(GemOm(ai,2))).*(2+(rand(size(GemOm(ai,2)))-0.5)/5),GemOm(ai,2),20,'r','filled')
% scatter(ones(size(GemOm(ai,3))).*(3+(rand(size(GemOm(ai,3)))-0.5)/5),GemOm(ai,3),20,'r','filled')
ylabel('Omissions (%)')
xlabel('SD (sec)');
ylim([0 70])

SPSS.SD.Om = [GemOm(ai,1), GemOm(ai,2), GemOm(ai,3)];
%GemOmStatGroup = [ones(size(GemOm(ai,1),1),1);  2*ones(size(GemOm(ai,2),1),1); 3*ones(size(GemOm(ai,3),1),1)];
%[OmStat{1},OmStat{2},OmStat{3}] = kruskalwallis(GemOmStat,GemOmStatGroup,'off');
%OmStatSD{4} = multcompare(OmStat{3},'Display','off');

subplot(2,4,3)
hold on
%notBoxPlot([GemPrem(ai,1) GemPrem(ai,2) GemPrem(ai,3)])
plot([(GemPrem(ai,3)./nrTrSD'*100)'; (GemPrem(ai,2)./nrTrSD'*100)'; (GemPrem(ai,1)./nrTrSD'*100)'],'k')
boxplot([GemPrem(ai,3)./nrTrSD'*100 GemPrem(ai,2)./nrTrSD'*100 GemPrem(ai,1)./nrTrSD'*100],'Labels',{'1.0','0.5','0.20'})
% scatter(ones(size(GemPrem(ai,1))).*(1+(rand(size(GemPrem(ai,1)))-0.5)/5),GemPrem(ai,1)./nrTrSD'*100,20,'r','filled')
% scatter(ones(size(GemPrem(ai,2))).*(2+(rand(size(GemPrem(ai,2)))-0.5)/5),GemPrem(ai,2)./nrTrSD'*100,20,'r','filled')
% scatter(ones(size(GemPrem(ai,3))).*(3+(rand(size(GemPrem(ai,3)))-0.5)/5),GemPrem(ai,3)./nrTrSD'*100,20,'r','filled')
ylabel('Premature responses (%)')
xlabel('SD (sec)');
ylim([0 20])
%ax = gca;
%ax.XAxis.TickLabel = {'0.2', '0.5', '1.0'};

SPSS.SD.Prem = [(GemPrem(ai,1)./nrTrSD'*100), (GemPrem(ai,2)./nrTrSD'*100), (GemPrem(ai,3)./nrTrSD'*100)];
%GemPremStatGroup = [ones(size(GemPrem(ai,1),1),1);  2*ones(size(GemPrem(ai,2),1),1); 3*ones(size(GemPrem(ai,3),1),1)];
%[PremStat{1},PremStat{2},PremStat{3}] = kruskalwallis(GemPremStat,GemPremStatGroup,'off');
%PremStatSD{4} = multcompare(PremStat{3},'Display','off');

subplot(2,4,4)
hold on
%notBoxPlot([GemPrem(ai,1) GemPrem(ai,2) GemPrem(ai,3)])
boxplot(nrTrSD)
scatter(ones(size(nrTrSD)).*(1+(rand(size(nrTrSD))-0.5)/5),nrTrSD,20,'k','filled')
%title('SD')
ylabel('Trials per session')
ylim([200 700])
xlabel('something')


%% Plot behavioral performance for VAR_ITI
for i = 1:size(Perf_ITI,1)
GemAcITI(i,1) = Perf_ITI{i,1}.accuracyPerITI(1,1);
GemAcITI(i,2) = Perf_ITI{i,1}.accuracyPerITI(1,2);
GemAcITI(i,3) = Perf_ITI{i,1}.accuracyPerITI(1,3);

GemOmITI(i,1) = Perf_ITI{i,1}.omissionPerITI(1,1);
GemOmITI(i,2) = Perf_ITI{i,1}.omissionPerITI(1,2);
GemOmITI(i,3) = Perf_ITI{i,1}.omissionPerITI(1,3);

GemPremITI(i,1) = Perf_ITI{i,1}.prematurePerITI(1,1);
GemPremITI(i,2) = Perf_ITI{i,1}.prematurePerITI(1,2);
GemPremITI(i,3) = Perf_ITI{i,1}.prematurePerITI(1,3);

GemCorITI(i,1) = numel(ITI.Correct{i, 1}.TrialPosITI50);
GemCorITI(i,2) = numel(ITI.Correct{i, 1}.TrialPosITI75);
GemCorITI(i,3) = numel(ITI.Correct{i, 1}.TrialPosITI125);

GemIncITI(i,1) = numel(ITI.Incorrect{i, 1}.TrialPosITI50);
GemIncITI(i,2) = numel(ITI.Incorrect{i, 1}.TrialPosITI75);
GemIncITI(i,3) = numel(ITI.Incorrect{i, 1}.TrialPosITI125);
end

% Find first position of unique pattern in variable ai
[~, ai, ~] = unique(GemAcITI(:,1),'stable');

% Make boxplots and do kruskall wallis statistics with multiple comparisons
subplot(2,4,5)
hold on
%notBoxPlot([GemAc(ai,1) GemAc(ai,2) GemAc(ai,3)])
plot([GemAcITI(ai,1)'; GemAcITI(ai,2)'; GemAcITI(ai,3)'],'k')
boxplot([GemAcITI(ai,1) GemAcITI(ai,2) GemAcITI(ai,3)],'Labels',{'5.0','7.5','12.5'})
% scatter(ones(size(GemAcITI(ai,1))).*(1+(rand(size(GemAcITI(ai,1)))-0.5)/5),GemAcITI(ai,1),20,'r','filled')
% scatter(ones(size(GemAcITI(ai,2))).*(2+(rand(size(GemAcITI(ai,2)))-0.5)/5),GemAcITI(ai,2),20,'r','filled')
% scatter(ones(size(GemAcITI(ai,3))).*(3+(rand(size(GemAcITI(ai,3)))-0.5)/5),GemAcITI(ai,3),20,'r','filled')
ylabel('Accuracy (%)')
xlabel('ITI (sec)');
ylim([50 100])

SPSS.ITI.Acc = [GemAcITI(ai,1), GemAcITI(ai,2), GemAcITI(ai,3)];
%GemAcStatGroup = [ones(size(GemAcITI(ai,1),1),1);  2*ones(size(GemAcITI(ai,2),1),1); 3*ones(size(GemAcITI(ai,3),1),1)];
%[AcStat{1},AcStat{2},AcStat{3}] = kruskalwallis(GemAcStat,GemAcStatGroup,'off');
%AcStatITI{4} = multcompare(AcStat{3},'Display','off');

subplot(2,4,6)
hold on
%notBoxPlot([GemOm(ai,1) GemOm(ai,2) GemOm(ai,3)])
plot([GemOmITI(ai,1)'; GemOmITI(ai,2)'; GemOmITI(ai,3)'],'k')
boxplot([GemOmITI(ai,1) GemOmITI(ai,2) GemOmITI(ai,3)],'Labels',{'5.0','7.5','12.5'})
% scatter(ones(size(GemOmITI(ai,1))).*(1+(rand(size(GemOmITI(ai,1)))-0.5)/5),GemOmITI(ai,1),20,'r','filled')
% scatter(ones(size(GemOmITI(ai,2))).*(2+(rand(size(GemOmITI(ai,2)))-0.5)/5),GemOmITI(ai,2),20,'r','filled')
% scatter(ones(size(GemOmITI(ai,3))).*(3+(rand(size(GemOmITI(ai,3)))-0.5)/5),GemOmITI(ai,3),20,'r','filled')
ylabel('Omissions (%)')
xlabel('ITI (sec)');
ylim([0 70])

for i = 1:size(nrTrialsStageITI,1)
    nrTrITI(i) = nrTrialsStageITI(i,1).TotalNr;
end

%%
SPSS.ITI.Prem = [(GemPremITI(ai,1)./nrTrITI'*100), (GemPremITI(ai,2)./nrTrITI'*100), (GemPremITI(ai,3)./nrTrITI'*100)];
%GemPremStatGroup = [ones(size(GemPremITI(ai,1),1),1);  2*ones(size(GemPremITI(ai,2),1),1); 3*ones(size(GemPremITI(ai,3),1),1)];
%[PremStat{1},PremStat{2},PremStat{3}] = kruskalwallis(GemPremStat,GemAcStatGroup,'off');
%PremStatITI{4} = multcompare(PremStat{3},'Display','off');

%%

SPSS.ITI.Om = [GemOmITI(ai,1), GemOmITI(ai,2), GemOmITI(ai,3)];
%GemOmStatGroup = [ones(size(GemOmITI(ai,1),1),1);  2*ones(size(GemOmITI(ai,2),1),1); 3*ones(size(GemOmITI(ai,3),1),1)];
%[OmStat{1},OmStat{2},OmStat{3}] = kruskalwallis(GemOmStat,GemOmStatGroup,'off');
%OmStatITI{4} = multcompare(OmStat{3},'Display','off');

subplot(2,4,7)
hold on
%notBoxPlot([GemPrem(ai,1) GemPrem(ai,2) GemPrem(ai,3)])
plot([(GemPremITI(ai,1)./nrTrITI'*100)'; (GemPremITI(ai,2)./nrTrITI'*100)'; (GemPremITI(ai,3)./nrTrITI'*100)'],'k')
boxplot([GemPremITI(ai,1)./nrTrITI'*100 GemPremITI(ai,2)./nrTrITI'*100 GemPremITI(ai,3)./nrTrITI'*100],'Labels',{'5.0','7.5','12.5'})
% scatter(ones(size(GemPremITI(ai,1))).*(1+(rand(size(GemPremITI(ai,1)))-0.5)/5),GemPremITI(ai,1)./nrTrITI'*100,20,'r','filled')
% scatter(ones(size(GemPremITI(ai,2))).*(2+(rand(size(GemPremITI(ai,2)))-0.5)/5),GemPremITI(ai,2)./nrTrITI'*100,20,'r','filled')
% scatter(ones(size(GemPremITI(ai,3))).*(3+(rand(size(GemPremITI(ai,3)))-0.5)/5),GemPremITI(ai,3)./nrTrITI'*100,20,'r','filled')
ylabel('Premature responses (%)')
xlabel('ITI (sec)');
ylim([0 20])
%ax = gca;
%ax.XAxis.TickLabel = {'5.0', '7.5', '12.5'};

subplot(2,4,8)
hold on
%notBoxPlot([GemPrem(ai,1) GemPrem(ai,2) GemPrem(ai,3)])
boxplot(nrTrITI)
scatter(ones(size(nrTrITI)).*(1+(rand(size(nrTrITI))-0.5)/5),nrTrITI,20,'k','filled')
%title('ITI')
ylabel('Trials per session')
ylim([200 700])
xlabel('something');

%% Make figure on Prem resp time for VAR_ITI
Lat.ITI50 = [];
Lat.ITI75 = [];
Lat.ITI125 = [];



n = 1;
for i = 1:size(PrematureITI,1)
    Lat.ITI50(n,:) = histcounts((PrematureITI(i).LatencyITI50),0:0.5:12.5);
    Lat.ITI75(n,:) = histcounts((PrematureITI(i).LatencyITI75),0:0.5:12.5);
    Lat.ITI125(n,:) = histcounts((PrematureITI(i).LatencyITI125),0:0.5:12.5);
    n = n+1;
end
for i = 1:size(Lat.ITI50,1)
    if Lat.ITI50(i,11) > 0
        Lat.ITI50(i,10) = Lat.ITI50(i,10)+Lat.ITI50(i,11);
        Lat.ITI50(i,11) = 0;
    elseif Lat.ITI75(i,16) > 0
        Lat.ITI75(i,15) = Lat.ITI75(i,15)+Lat.ITI50(i,16);
        Lat.ITI75(i,16) = 0;
    end
end
x_ITI125 = 0.5:0.5:12.5;

figure
bar(x_ITI125,[mean(Lat.ITI50,1)' mean(Lat.ITI75,1)' mean(Lat.ITI125,1)'])
title('Premature response time Var. ITI')
xlabel('Time (sec)')
ylabel('Responses/session (mean)')
legend({'ITI 5.0', 'ITI 7.5', 'ITI 12.5'})

% Find neuron indexes to use to find cue orientation moment
ITIIx = [];
for l = 1:numel(VAR_ITIix)
    for i = 1:size(All{5, 1},1)
        if All{5, 1}{i,2} == MED_PC{VAR_ITIix(l),1} & contains(All{5, 1}{i,1}, 'ITI')
            ITIIx(end+1) = i;
            break
        end
    end
end






% for i = 1:size(PrematureITI,1)
%     PreLatITI_125(i,:) = histcounts(PrematureITI(i).LatencyITI125, 0:0.5:12.5, 'Normalization', 'probability');
%     PreLatITIgood_125(i,:) = histcounts(PrematureITI(i).LatencyITI125(VideoTrialsRemoval_ITI125{i,1}.Prem_ThreeSecsIx), 0:0.5:12.5, 'Normalization', 'probability');
%     PreLatITI_75(i,:) = histcounts(PrematureITI(i).LatencyITI75, 0:0.5:7.5, 'Normalization', 'probability');
%     PreLatITIgood_75(i,:) = histcounts(PrematureITI(i).LatencyITI75(VideoTrialsRemoval_ITI75{i,1}.Prem_ThreeSecsIx), 0:0.5:7.5, 'Normalization', 'probability');
%     PreLatITI_50(i,:) = histcounts(PrematureITI(i).LatencyITI50, 0:0.5:5.0, 'Normalization', 'probability');
%     PreLatITIgood_50(i,:) = histcounts(PrematureITI(i).LatencyITI50(VideoTrialsRemoval_ITI50{i,1}.Prem_ThreeSecsIx), 0:0.5:5.0, 'Normalization', 'probability');
%     CorTeshold_125(i,:) = histcounts(eventVideoBehavCor_ITI125{i, 1}.treshold', 0:0.5:12.5, 'Normalization', 'probability');
%     PremTeshold_125(i,:) = histcounts(eventVideoBehavPrem_ITI125{i, 1}.treshold'  , 0:0.5:12.5, 'Normalization', 'probability');
%         CorTeshold_75(i,:) = histcounts(eventVideoBehavCor_ITI75{i, 1}.treshold', 0:0.5:7.5, 'Normalization', 'probability');
%     PremTeshold_75(i,:) = histcounts(eventVideoBehavPrem_ITI75{i, 1}.treshold'  , 0:0.5:7.5, 'Normalization', 'probability');
%         CorTeshold_50(i,:) = histcounts(eventVideoBehavCor_ITI50{i, 1}.treshold', 0:0.5:5.0, 'Normalization', 'probability');
%     try
%         PremTeshold_50(i,:) = histcounts(eventVideoBehavPrem_ITI50{i, 1}.treshold'  , 0:0.5:5.0, 'Normalization', 'probability');
%     catch
%         PremTeshold_50(i,:) = nan;
%     end
% end

% for i = 1:size(PrematureSD,1)
%     %PreLatSD(i,:) = histcounts(PrematureSD(i).LatencyITI5, 0:0.5:5.0, 'Normalization', 'probability');
%     %PreLatSDgood(i,:) = histcounts(PrematureITI(i).LatencyITI5(VideoTrialsRemoval_SD{i,1}.Cor_ThreeSecsIx), 0:0.5:5.0, 'Normalization', 'probability');
%     Treshold.CorTresholdSD(i,:) = histcounts(eventVideoBehavCor_SD{i, 1}.treshold , 0:0.5:5.0, 'Normalization', 'count'); 
%     Treshold.IncTresholdSD(i,:) = histcounts(eventVideoBehavInc_SD{i, 1}.treshold , 0:0.5:5.0, 'Normalization', 'count'); 
%     Treshold.OmTresholdSD(i,:) = histcounts(eventVideoBehavOm_SD{i, 1}.treshold , 0:0.5:5.0, 'Normalization', 'count'); 
%     Treshold.PremTresholdSD(i,:) = histcounts(eventVideoBehavPrem_SD{i, 1}.treshold , 0:0.5:5.0, 'Normalization', 'count');
%     Treshold.CorTresholdSD(i,:) = Treshold.CorTresholdSD(i,:)./sum(Treshold.CorTresholdSD(i,:));
%     Treshold.IncTresholdSD(i,:) = Treshold.IncTresholdSD(i,:)./sum(Treshold.IncTresholdSD(i,:));
%     Treshold.OmTresholdSD(i,:) = Treshold.OmTresholdSD(i,:)./sum(Treshold.OmTresholdSD(i,:));
%     Treshold.PremTresholdSD(i,:) = Treshold.PremTresholdSD(i,:)./sum(Treshold.PremTresholdSD(i,:));
% end

% figure
% hold on
% for pos = 1:size(Treshold.CorTresholdSD,1)
% plot(0.5:0.5:5,Treshold.CorTresholdSD(pos,:),'k')
% end
% for pos = 1:size(Treshold.IncTresholdSD,1)
% plot(0.5:0.5:5,Treshold.IncTresholdSD(pos,:),'r')
% end
% for pos = 1:size(Treshold.OmTresholdSD,1)
% plot(0.5:0.5:5,Treshold.OmTresholdSD(pos,:),'b')
% end
% plot(0.5:0.5:5,mean(Treshold.CorTresholdSD),'Color','k','LineWidth',3)
% plot(0.5:0.5:5,mean(Treshold.IncTresholdSD),'Color','r','LineWidth',3)
% plot(0.5:0.5:5,mean(Treshold.OmTresholdSD),'Color','b','LineWidth',3)

% figure
% hold on
% SEM = nanstd(Treshold.CorTresholdSD,1);%/sqrt(size(psthBinsValue.Cor.trial_start,1)); 
% errorbar(0.5:0.5:5,mean(Treshold.CorTresholdSD),SEM,'Color','k')
% SEM = nanstd(Treshold.IncTresholdSD,1);%/sqrt(size(psthBinsValue.Cor.trial_start,1)); 
% errorbar(0.5:0.5:5,mean(Treshold.IncTresholdSD),SEM,'Color','r')
% SEM = nanstd(Treshold.OmTresholdSD,1);%/sqrt(size(psthBinsValue.Cor.trial_start,1)); 
% errorbar(0.5:0.5:5,mean(Treshold.OmTresholdSD),SEM,'Color','b')
% SEM = nanstd(Treshold.PremTresholdSD,1);%/sqrt(size(psthBinsValue.Cor.trial_start,1)); 
% errorbar(0.5:0.5:5,mean(Treshold.PremTresholdSD),SEM,'Color','c')
% xlabel('Timebin (sec)')
% ylabel('Probability')
% title('Start trial to cue orientation time VAR_SD')
% legend({'Correct','Incorrect','Omission'})
% 
% figure
% hold on
% SEM = nanstd(Treshold.CorTresholdITI,1);%/sqrt(size(psthBinsValue.Cor.trial_start,1)); 
% errorbar(0.5:0.5:5,mean(Treshold.CorTresholdITI),SEM,'Color','k')
% SEM = nanstd(Treshold.IncTresholdITI,1);%/sqrt(size(psthBinsValue.Cor.trial_start,1)); 
% errorbar(0.5:0.5:5,mean(Treshold.IncTresholdITI),SEM,'Color','r')
% SEM = nanstd(Treshold.OmTresholdITI,1);%/sqrt(size(psthBinsValue.Cor.trial_start,1)); 
% errorbar(0.5:0.5:5,mean(Treshold.OmTresholdITI),SEM,'Color','b')
% SEM = nanstd(Treshold.PremTresholdITI,1);%/sqrt(size(psthBinsValue.Cor.trial_start,1)); 
% errorbar(0.5:0.5:5,mean(Treshold.PremTresholdITI),SEM,'Color','c')
% xlabel('Timebin (sec)')
% ylabel('Probability')
% title('Start trial to cue orientation time VAR_ITI')
% legend({'Correct','Incorrect','Omission'})
% 
% figure
% subplot(1,3,1)
% plot(mean(PreLatITI_125))
% hold on
% plot(nanmean(PreLatITIgood_125))
% subplot(1,3,2)
% plot(mean(PreLatITI_75))
% hold on
% plot(nanmean(PreLatITIgood_75))
% subplot(1,3,3)
% plot(nanmean(PreLatITI_50))
% hold on
% plot(nanmean(PreLatITIgood_50))
