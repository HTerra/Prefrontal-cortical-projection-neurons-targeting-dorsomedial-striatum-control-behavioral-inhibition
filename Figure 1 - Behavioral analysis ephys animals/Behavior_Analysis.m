function [SD, ITI, SPSS, nrTrialsStageITI, nrTrialsStageSD, PrematureITI, PrematureSD, MED_PC] = Behavior_Analysis(All,neuronIndex)

% Plots Figure 1

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

SPSS.SD.Acc = [GemAc(ai,1), GemAc(ai,2), GemAc(ai,3)];

subplot(2,4,2)
hold on
plot([GemOm(ai,3)'; GemOm(ai,2)'; GemOm(ai,1)'],'k')
boxplot([GemOm(ai,3) GemOm(ai,2) GemOm(ai,1)],'Labels',{'1.0','0.5','0.20'})
ylabel('Omissions (%)')
xlabel('SD (sec)');
ylim([0 70])

SPSS.SD.Om = [GemOm(ai,1), GemOm(ai,2), GemOm(ai,3)];

subplot(2,4,3)
hold on
plot([(GemPrem(ai,3)./nrTrSD'*100)'; (GemPrem(ai,2)./nrTrSD'*100)'; (GemPrem(ai,1)./nrTrSD'*100)'],'k')
boxplot([GemPrem(ai,3)./nrTrSD'*100 GemPrem(ai,2)./nrTrSD'*100 GemPrem(ai,1)./nrTrSD'*100],'Labels',{'1.0','0.5','0.20'})
ylabel('Premature responses (%)')
xlabel('SD (sec)');
ylim([0 20])

SPSS.SD.Prem = [(GemPrem(ai,1)./nrTrSD'*100), (GemPrem(ai,2)./nrTrSD'*100), (GemPrem(ai,3)./nrTrSD'*100)];

subplot(2,4,4)
hold on
boxplot(nrTrSD)
scatter(ones(size(nrTrSD)).*(1+(rand(size(nrTrSD))-0.5)/5),nrTrSD,20,'k','filled')
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
plot([GemAcITI(ai,1)'; GemAcITI(ai,2)'; GemAcITI(ai,3)'],'k')
boxplot([GemAcITI(ai,1) GemAcITI(ai,2) GemAcITI(ai,3)],'Labels',{'5.0','7.5','12.5'})
ylabel('Accuracy (%)')
xlabel('ITI (sec)');
ylim([50 100])

SPSS.ITI.Acc = [GemAcITI(ai,1), GemAcITI(ai,2), GemAcITI(ai,3)];

subplot(2,4,6)
hold on
plot([GemOmITI(ai,1)'; GemOmITI(ai,2)'; GemOmITI(ai,3)'],'k')
boxplot([GemOmITI(ai,1) GemOmITI(ai,2) GemOmITI(ai,3)],'Labels',{'5.0','7.5','12.5'})
ylabel('Omissions (%)')
xlabel('ITI (sec)');
ylim([0 70])

for i = 1:size(nrTrialsStageITI,1)
    nrTrITI(i) = nrTrialsStageITI(i,1).TotalNr;
end

%%
SPSS.ITI.Prem = [(GemPremITI(ai,1)./nrTrITI'*100), (GemPremITI(ai,2)./nrTrITI'*100), (GemPremITI(ai,3)./nrTrITI'*100)];

%%

SPSS.ITI.Om = [GemOmITI(ai,1), GemOmITI(ai,2), GemOmITI(ai,3)];

subplot(2,4,7)
hold on
plot([(GemPremITI(ai,1)./nrTrITI'*100)'; (GemPremITI(ai,2)./nrTrITI'*100)'; (GemPremITI(ai,3)./nrTrITI'*100)'],'k')
boxplot([GemPremITI(ai,1)./nrTrITI'*100 GemPremITI(ai,2)./nrTrITI'*100 GemPremITI(ai,3)./nrTrITI'*100],'Labels',{'5.0','7.5','12.5'})
ylabel('Premature responses (%)')
xlabel('ITI (sec)');
ylim([0 20])

subplot(2,4,8)
hold on
boxplot(nrTrITI)
scatter(ones(size(nrTrITI)).*(1+(rand(size(nrTrITI))-0.5)/5),nrTrITI,20,'k','filled')
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

SDIx = [];
for l = 1:numel(VAR_SDix)
    for i = 1:size(All{5, 1},1)
        if All{5, 1}{i,2} == MED_PC{VAR_SDix(l),1} & contains(All{5, 1}{i,1}, 'SD')
            SDIx(end+1) = i;
            break
        end
    end
end

%% performance index per animal
ITIixanimal = {};
for i = 1:size(ITIIx,2)
ITIixanimal{i} = All{8, 1}{ITIIx(i),1};
end
animalIxes = {1:5, 6:13, 14:19, 20:24, 25:30};
SPSS.ITI.GemAc_Animal = [];
SPSS.ITI.GemOm_Animal = [];
SPSS.ITI.GemPrem_Animal = [];
GemAcITI_Animal = [];
GemPremITI_Animal = [];
GemOmITI_Animal = [];
for i = 1:5
SPSS.ITI.GemAc_Animal(i,:) = mean(SPSS.ITI.Acc(animalIxes{i},:),1);
SPSS.ITI.GemPrem_Animal(i,:) = mean(SPSS.ITI.Prem(animalIxes{i},:),1);
SPSS.ITI.GemOm_Animal(i,:) = mean(SPSS.ITI.Om(animalIxes{i},:),1);
GemAcITI_Animal(i,:) = mean(SPSS.ITI.Acc(animalIxes{i},:),1);
GemPremITI_Animal(i,:) = mean(SPSS.ITI.Prem(animalIxes{i},:),1);
GemOmITI_Animal(i,:) = mean(SPSS.ITI.Om(animalIxes{i},:),1);
end

SDixanimal = {};
for i = 1:size(SDIx,2)
SDixanimal{i} = All{8, 1}{SDIx(i),1};
end
animalIxes = {1:7, 8:15, 16:20, 21:25, 26:33};
SPSS.SD.GemAc_Animal = [];
SPSS.SD.GemOm_Animal = [];
SPSS.SD.GemPrem_Animal = [];
GemAc_Animal = [];
GemOm_Animal = [];
GemPrem_Animal = [];
SPSS.SD.GemAc_Animal = [];
SPSS.SD.GemPrem_Animal = [];
SPSS.SD.GemOm_Animal = [];
for i = 1:5
SPSS.SD.GemAc_Animal(i,:) = mean(SPSS.SD.Acc(animalIxes{i},:),1);
SPSS.SD.GemPrem_Animal(i,:) = mean(SPSS.SD.Prem(animalIxes{i},:),1);
SPSS.SD.GemOm_Animal(i,:) = mean(SPSS.SD.Om(animalIxes{i},:),1);
GemAc_Animal(i,:) = mean(SPSS.SD.Acc(animalIxes{i},:),1);
GemPrem_Animal(i,:) = mean(SPSS.SD.Prem(animalIxes{i},:),1);
GemOm_Animal(i,:) = mean(SPSS.SD.Om(animalIxes{i},:),1);
end
%

figure
subplot(2,4,1)
hold on
plot([GemAc_Animal(:,3)'; GemAc_Animal(:,2)'; GemAc_Animal(:,1)'],'k')
boxplot([GemAc_Animal(:,3) GemAc_Animal(:,2) GemAc_Animal(:,1)],'Labels',{'1.0','0.5','0.20'})
ylabel('Accuracy (%)')
xlabel('SD (sec)');
ylim([50 100])

subplot(2,4,2)
hold on
plot([GemOm_Animal(:,3)'; GemOm_Animal(:,2)'; GemOm_Animal(:,1)'],'k')
boxplot([GemOm_Animal(:,3) GemOm_Animal(:,2) GemOm_Animal(:,1)],'Labels',{'1.0','0.5','0.20'})
ylabel('Omissions (%)')
xlabel('SD (sec)');
ylim([0 70])

subplot(2,4,3)
hold on
plot([SPSS.SD.GemPrem_Animal(:,3)'; SPSS.SD.GemPrem_Animal(:,2)'; SPSS.SD.GemPrem_Animal(:,1)'],'k')
boxplot([SPSS.SD.GemPrem_Animal(:,3)'; SPSS.SD.GemPrem_Animal(:,2)'; SPSS.SD.GemPrem_Animal(:,1)']','Labels',{'1.0','0.5','0.20'})
ylabel('Premature responses (%)')
xlabel('SD (sec)');
ylim([0 20])

subplot(2,4,5)
hold on
plot([GemAcITI_Animal(:,1)'; GemAcITI_Animal(:,2)'; GemAcITI_Animal(:,3)'],'k')
boxplot([GemAcITI_Animal(:,1) GemAcITI_Animal(:,2) GemAcITI_Animal(:,3)],'Labels',{'5.0','7.5','12.5'})
ylabel('Accuracy (%)')
xlabel('ITI (sec)');
ylim([50 100])

subplot(2,4,6)
hold on
plot([GemOmITI_Animal(:,1)'; GemOmITI_Animal(:,2)'; GemOmITI_Animal(:,3)'],'k')
boxplot([GemOmITI_Animal(:,1) GemOmITI_Animal(:,2) GemOmITI_Animal(:,3)],'Labels',{'5.0','7.5','12.5'})
ylabel('Omissions (%)')
xlabel('ITI (sec)');
ylim([0 70])

subplot(2,4,7)
hold on
plot([SPSS.ITI.GemPrem_Animal(:,1)'; SPSS.ITI.GemPrem_Animal(:,2)'; SPSS.ITI.GemPrem_Animal(:,3)'],'k')
boxplot([SPSS.ITI.GemPrem_Animal(:,1)'; SPSS.ITI.GemPrem_Animal(:,2)'; SPSS.ITI.GemPrem_Animal(:,3)']','Labels',{'5.0','7.5','12.5'})
ylabel('Premature responses (%)')
xlabel('ITI (sec)');
ylim([0 20])
