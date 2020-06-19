function [Correct, Incorrect, Omission, Premature, Perserverative, nrTrialsStage, Perf]  = CC_ephys_var_ITI_SD_TTL_Master(stage,timeBlocks,MED_PC,varargin)

% Function used to get #th of trial per ITI or SD condition per type of response.
% Is used to split TTL pulses in Neurorighter per ITI or SD condition
% (with get_eventTimesNeurorighter.m). Can also be used to plot performance
% per animal.

% Input:
% - stage = stage to analyse with var ITI or SD
% - timeBlocks = Number of blocks the 2.5 hour session should be divided in.
%  - varargin:
%     'plotTimeBlocks' plots the number of trials and performance per defined timeblock (timeBlocks) in a bar graph.
%     'plotVarITI' plots the performance per ITI.
%     'plotVarSD' plots the performance per SD.


% Ouput:
% Correct   = Struct containing: Response indexes divided by ITI and SD condition,
%           Total # of trials Correct, Position of correct trials for all
%           trials done.
% Incorrect = Struct containing: Response indexes divided by ITI and SD condition,
%           Total # of trials Correct, Position of correct trials for all
%           trials done.
% Omission  = Struct containing: Response indexes divided by ITI and SD condition,
%           Total # of trials Correct, Position of correct trials for all
%           trials done.
% Premature = Struct containing: Response indexes divided by ITI and SD condition,
%           Total # of trials Correct, Position of correct trials for all
%           trials done.
% Different plots

% V1.0
% h.terra@vu.nl
% 25th of June 2017

%% Get MED-PC output file in .txt form.
%RecParam;

%% Find where x array starts and imports the whole X array.
[data,delimiterOut,headerlinesOut]=importdata(MED_PC,'',60);
idx = (find(ismember(data(:,1), 'X:')))+60-headerlinesOut;
clear data delimiterOut headerlinesOut
data=importdata(MED_PC,' ',idx);

%% Cut the X array off and the last trial.
endFind = 1: 5: length(data.data);

for basta = endFind
    if data.data(basta,1) == 0
        final = basta;
        data.data = data.data(1:final-1,:);
        break
    end
end

data = data.data;

clearvars -except data MED_PC stage timeBlocks varargin

%% Divide in data matrix per stage.
indexing = 1:5:length(data)-9;
stageUpIndex = 1;
n = 2;
for i = indexing
    if data(i+5,3) - data(i,3) > 0
        stageUpIndex(n) = i+5;
        n = n+1;
    end
end

clearvars -except data MED_PC stageUpIndex stage timeBlocks varargin

n = 1;
for stages = stageUpIndex
    temp = sprintf('stage%i', data(stages,3));
    if n == length(stageUpIndex)
        temp2 = length(data);
    else
        temp2 = stageUpIndex(n+1)-1;
    end
    data2.(temp) = data(stages:temp2,:);
        if data(stages,3) == stage
            day = n;
        end
    n = n+1;
end

clearvars -except data data2 MED_PC stageUpIndex stage timeBlocks varargin day

stageName = sprintf('stage%i', stage);
try
    data = data2.(stageName);
catch
    warning('Stage not present. Please indicate other stage')
end

%Normalize all times to the start of the first trial that day.
%if day > 1
    startTrial = 1:5:length(data);
    respIndex = 4:5:length(data);
    
    data(respIndex,2) = data(respIndex,2) - data(1,1);
    data(startTrial,1) = data(startTrial,1) - data(1,1);
%end

clear data2 temp startTrial respIndex

%% Make block analysis. Blocksize indicated by "timeBLocks".

pos = 1:5:length(data);
blockSize = (150*60)/timeBlocks;
blockStartIndex = ones(1,timeBlocks);

for i = 1:timeBlocks
    temp = sprintf('timeBlock%i', i);
    dataBlock.(temp) = [];
end

n = 1;
for i = pos;
    temp = sprintf('timeBlock%i', n);
    if data(i,1) > n*blockSize
        blockStartIndex(n+1) = i;
        dataBlock.(temp) = data(blockStartIndex(n):blockStartIndex(n+1)-1,:);
        n = n+1;
    else n == timeBlocks;
        dataBlock.(temp) = data(blockStartIndex(n):end,:);
    end
end

%% Give struct with values that indicate the index of different ITI and SD durations per trial type (correct, incorrect, omission, premature, perserverative). Used for linking neural activity per ITI and SD condition.

[Correct, Incorrect, Omission, Premature, Perserverative, nrTrialsStage] = CC_ephys_var_ITI_SD_TTL_Child(data);

% Index per timeblock.
try
    for block = 1:timeBlocks
        temp = sprintf('timeBlock%i', block);
        [Correct.(temp), Incorrect.(temp), Omission.(temp), Premature.(temp), Perserverative.(temp), nrTrialsStage.(temp)] = CC_ephys_var_ITI_SD_TTL_Child(dataBlock.(temp),block);
    end
catch
    formatSpec = 'timeblock %d %s';
    A1 = 'is empty. Consider smaller size of timeBlocks';
    warning(formatSpec,block,A1);
end


%% Calc performance for VAR-ITI sessions
Perf.performancePerITI = zeros(3,3);
Perf.accuracyPerITI = zeros(1,3);
Perf.omissionPerITI = zeros(1,3);
Perf.prematurePerITI = zeros(1,3);
Perf.perserverativePerITI = zeros(1,3);

pos = 2:5:length(data);
ITI = unique(data(pos,1));
for l = 1:length(ITI)
    s = num2str(ITI(l)*10);
    ITIname = sprintf('TrialPosITI%s', s);
    n = 1;
    Perf.performancePerITI(l,n) = length(Correct.(ITIname));
    n = n+1;
    Perf.performancePerITI(l,n) = length(Incorrect.(ITIname));
    n = n+1;
    Perf.performancePerITI(l,n) = length(Omission.(ITIname));
    Perf.accuracyPerITI(l) = (Perf.performancePerITI(l,n-2) / (Perf.performancePerITI(l,n-2) + Perf.performancePerITI(l,n-1)))*100;
    Perf.omissionPerITI(l) = (Perf.performancePerITI(l,n) / (Perf.performancePerITI(l,n-2) + Perf.performancePerITI(l,n-1) + Perf.performancePerITI(l,n)))*100;
    Perf.prematurePerITI(l) = length(Premature.(ITIname));
    Perf.perserverativePerITI(l) = length(Perserverative.(ITIname));
end

%% Calc performance for VAR-SD sessions
Perf.performancePerSD = zeros(3,3);
Perf.accuracyPerSD = zeros(1,3);
Perf.omissionPerSD = zeros(1,3);
Perf.prematurePerSD = zeros(1,3);
Perf.perserverativePerSD = zeros(1,3);

pos = 2:5:length(data);
SD = unique(data(pos,5));
for l = 1:length(SD)
    s = num2str(SD(l)*10);
    SDname = sprintf('TrialPosSD%s', s);
    n = 1;
    Perf.performancePerSD(l,n) = length(Correct.(SDname));
    n = n+1;
    Perf.performancePerSD(l,n) = length(Incorrect.(SDname));
    n = n+1;
    Perf.performancePerSD(l,n) = length(Omission.(SDname));
    Perf.accuracyPerSD(l) = (Perf.performancePerSD(l,n-2) / (Perf.performancePerSD(l,n-2) + Perf.performancePerSD(l,n-1)))*100;
    Perf.omissionPerSD(l) = (Perf.performancePerSD(l,n) / (Perf.performancePerSD(l,n-2) + Perf.performancePerSD(l,n-1) + Perf.performancePerSD(l,n)))*100;
    Perf.prematurePerSD(l) = length(Premature.(SDname));
    Perf.perserverativePerSD(l) = length(Perserverative.(SDname));
end

%% add plotting options

if nargin;
    for iEvent=1:size(varargin,2)
        switch varargin{iEvent}
            case 'plotTimeBlocks'
                performancePerBlock = zeros(timeBlocks,3);
                accuracyPerBlock = zeros(1,timeBlocks);
                omissionPerBlock = zeros(1,timeBlocks);

                for block = 1:timeBlocks;
                    n = 1;
                    temp = sprintf('timeBlock%i', block);
                    trialsPerBlock(block) = nrTrialsStage.(temp).TotalNr;
                    performancePerBlock(block,n) = Correct.(temp).TotalNr;
                    n = n+1;
                    performancePerBlock(block,n) = Incorrect.(temp).TotalNr;
                    n = n+1;
                    performancePerBlock(block,n) = Omission.(temp).TotalNr;
                    accuracyPerBlock(block) = (performancePerBlock(block,n-2) / (performancePerBlock(block,n-2) + performancePerBlock(block,n-1)))*100;
                    omissionPerBlock(block) = (performancePerBlock(block,n) / (performancePerBlock(block,n-2) + performancePerBlock(block,n-1) + performancePerBlock(block,n)))*100;
                end
                
                figure
                subplot(2,2,1);
                bar(accuracyPerBlock)
                title('Accuracy per timeblock')
                ylabel('Accuracy (%)');

                subplot(2,2,2);
                bar(omissionPerBlock)
                title('Omissions per timeblock')
                ylabel('Omissions (%)');

                subplot(2,2,[3, 4]);
                bar(performancePerBlock, 'stacked')
                title('Overall performance per timeblock')
                xlabel('Timeblock');
                ylabel('Trials');
                legend('Correct', 'Incorrect', 'Omissions');

                figure
                bar(trialsPerBlock)
                title('Number of trials per timeblock');
                xlabel('Block');
                ylabel('Trials');
                
            case 'plotVarITI'

                pos = 2:5:length(data);
                ITI = unique(data(pos,1));
                
                for i = 1:length(ITI);
                    s = num2str(ITI(i));
                    ITIlabels{i} = s;
                end
                
                figure
                subplot(2,2,1);
                bar(Perf.accuracyPerITI)
                title('Accuracy per ITI')
                ylabel('Accuracy (%)');
                set(gca,'XTickLabel', ITIlabels)

                subplot(2,2,2);
                bar(Perf.omissionPerITI)
                title('Omissions per ITI')
                ylabel('Omissions (%)');
                set(gca,'XTickLabel', ITIlabels)

                subplot(2,2,[3, 4]);
                bar(Perf.performancePerITI, 'stacked')
                title('Overall performance per ITI')
                xlabel('Inter-trial interval (sec)');
                ylabel('Trials');
                legend('Correct', 'Incorrect', 'Omissions');
                set(gca,'XTickLabel', ITIlabels)
                
                figure
                subplot(2,2,1);
                bar(Perf.prematurePerITI)
                title('Number of premature trials per ITI');
                xlabel('Inter-trial interval (sec)');
                ylabel('Trials');
                set(gca,'XTickLabel', ITIlabels)
                
                subplot(2,2,2);
                bar(Perf.perserverativePerITI)
                title('Number of Perserverative trials per ITI')
                xlabel('Inter-trial interval (sec)');
                ylabel('Trials');
                set(gca,'XTickLabel', ITIlabels)
                
            case 'plotVarSD'
                performancePerSD = zeros(3,3);
                accuracyPerSD = zeros(1,3);
                omissionPerSD = zeros(1,3);
                prematurePerSD = zeros(1,3);
                perserverativePerSD = zeros(1,3);

                pos = 2:5:length(data);
                SD = unique(data(pos,5));
                for l = 1:length(SD)
                    s = num2str(SD(l)*10);
                    SDname = sprintf('TrialPosSD%s', s);
                    n = 1;
                    performancePerSD(l,n) = length(Correct.(SDname));
                    n = n+1;
                    performancePerSD(l,n) = length(Incorrect.(SDname));
                    n = n+1;
                    performancePerSD(l,n) = length(Omission.(SDname));
                    accuracyPerSD(l) = (performancePerSD(l,n-2) / (performancePerSD(l,n-2) + performancePerSD(l,n-1)))*100;
                    omissionPerSD(l) = (performancePerSD(l,n) / (performancePerSD(l,n-2) + performancePerSD(l,n-1) + performancePerSD(l,n)))*100;
                    prematurePerSD(l) = length(Premature.(SDname));
                    perserverativePerSD(l) = length(Perserverative.(SDname));
                end
                
                for i = 1:length(SD);
                    s = num2str(SD(i));
                    SDlabels{i} = s;
                end
                
                figure
                subplot(2,2,1);
                bar(accuracyPerSD)
                title('Accuracy per SD')
                ylabel('Accuracy (%)');
                set(gca,'XTickLabel', SDlabels)

                subplot(2,2,2);
                bar(omissionPerSD)
                title('Omissions per SD');
                ylabel('Omissions (%)');
                set(gca,'XTickLabel', SDlabels)

                subplot(2,2,[3, 4]);
                bar(performancePerSD, 'stacked')
                title('Overall performance per SD')
                xlabel('Stimulus duration (sec)');
                ylabel('Trials');
                legend('Correct', 'Incorrect', 'Omissions');
                set(gca,'XTickLabel', SDlabels)
                
                figure
                subplot(2,2,1);
                bar(prematurePerSD)
                title('Number of premature trials per SD');
                xlabel('Stimulus duration (sec)');
                ylabel('Trials');
                set(gca,'XTickLabel', SDlabels)
                
                subplot(2,2,2);
                bar(perserverativePerSD)
                title('Number of Perserverative trials per SD')
                xlabel('Stimulus duration (sec)');
                ylabel('Trials');
                set(gca,'XTickLabel', SDlabels)
            otherwise
                warning('Unknown input argument')
        end
    end
end
