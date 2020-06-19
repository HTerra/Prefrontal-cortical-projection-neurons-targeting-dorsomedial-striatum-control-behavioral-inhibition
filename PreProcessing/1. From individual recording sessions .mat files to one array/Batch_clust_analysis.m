% Scripts does analysis on batch cluster data. Requires that all cluster
% preprocessing is done.

%% Load in the seperate session files and place in data array "All"
files = dir('*.mat');
n = 1;
for i = 1:length(files)
    file{n,1} = files(i,1).name;
    n = n+1;
end

clear files n

load(file{1});

% store field names
names = fieldnames(cluster);

All = struct2cell(cluster);
for veld = [16 17 18];
    All{veld,1}{1,end+11-size(All{veld,1},2)} = {};
end
for veld = 19;
        All{veld,1} = {};
end
for l = 1:size(All{1,1},1)
    for i = 1:size(All{1,1}{l,1},1)
        All{23,1}{l,1}(i,:) = cat(2,All{1,1}{l,1}(i,:),All{18,1}{l,i});
    end
end

All{1,1} = All{23,1};

clear cluster

% load all struct and make in array with struct2cell

for f = 2:length(file)
    load(file{f});
    clusterCell = struct2cell(cluster);
    for veld = [16 17 18];
        clusterCell{veld,1}{1,end+11-size(clusterCell{veld,1},2)} = {};
    end
    
    % Remove waveforms
    for veld = 19;
        clusterCell{veld,1} = {};
    end

    for l = 1:size(clusterCell{1,1},1)
        for i = 1:size(clusterCell{1,1}{l,1},1)
            clusterCell{23,1}{l,1}(i,:) = cat(2,clusterCell{1,1}{l,1}(i,:),clusterCell{18,1}{l,i});
        end
    end

    clusterCell{1,1} = clusterCell{23,1};
    
    for fi = 1:length(All)
        try
            All{fi,1} = cat(1, All{fi,1}, clusterCell{fi,1});
        catch
            All{fi,1} = nan;
        end
    end
    clear cluster
end

All = All(1:22,1);
All(:,2) = names;

clearvars -except All

ClustIx = 1:size(All{1,1},1);

%% calculate waveform amplitude during behavior and store at All{15,1}(:,9) (9th position in qualityWVBehav)
All{3,1}(:,8) = nan;

for i = 1:size(ClustIx,2)
    All{3,1}(ClustIx(i),8) = abs(min(All{13,1}(ClustIx(i),:)));
end

%% Add timestamps of Cue on with trial type for Correct, Incorrect, Omission, Prem and Pers in CH17, CH18, CH19, Ch20 and CH21. Add timestamps of Start trial
% for NeuPos = 1:size(All{22,1},1)
%     CueOn = [All{22,1}{NeuPos,1}.CH1 All{22,1}{NeuPos,1}.CH2 All{22,1}{NeuPos,1}.CH3 All{22,1}{NeuPos,1}.CH4 All{22,1}{NeuPos,1}.CH5];
%     for i = 1:size(All{22,1}{NeuPos,1}.CH6,2)
%         [cCor(i) indexCor(i)] = min(abs(CueOn-All{22,1}{NeuPos,1}.CH6(i)));
%     end
%     for i = 1:size(All{22,1}{NeuPos,1}.CH8,2)
%         [cInc(i) indexInc(i)] = min(abs(CueOn-All{22,1}{NeuPos,1}.CH8(i)));
%     end
%     for i = 1:size(All{22,1}{NeuPos,1}.CH10,2)
%         [cOm(i) indexOm(i)] = min(abs(CueOn-All{22,1}{NeuPos,1}.CH10(i)));
%     end
%     for i = 1:size(All{22,1}{NeuPos,1}.CH12,2)
%         [cPrem(i) indexPrem(i)] = min(abs(CueOn-All{22,1}{NeuPos,1}.CH12(i)));
%     end
%     for i = 1:size(All{22,1}{NeuPos,1}.CH13,2)
%         [cPers(i) indexPers(i)] = min(abs(CueOn-All{22,1}{NeuPos,1}.CH13(i)));
%     end
%     All{22,1}{NeuPos,1}.CH17 = All{22,1}{NeuPos,1}.CH6-cCor;
%     All{22,1}{NeuPos,1}.CH18 = All{22,1}{NeuPos,1}.CH8-cInc;
%     All{22,1}{NeuPos,1}.CH19 = All{22,1}{NeuPos,1}.CH10-cOm;
%     All{22,1}{NeuPos,1}.CH20 = All{22,1}{NeuPos,1}.CH12-cPrem;
%     All{22,1}{NeuPos,1}.CH21 = All{22,1}{NeuPos,1}.CH13-cPers;
%     clear cCor indexCor cInc indexInc cOm indexOm cPrem indexPrem cPers indexPers
% end

% Insert extra behavioral event cues like, cue on, mag etc. Also add extra
% behav info like latencies in the eventTimes section.
All = insertBehavVars(All);

%% Redo SALT analysis (including latency,jitter, reliability and waveformcompare)
for neuron = 1:size(All{1, 1},1)
    for stimBlock = 1:size(All{1, 1}{neuron, 1},1)
        [spt_test, spt_baseline, FSLatency, jitter, reliability, spt_spikeIx, baselineTimeArray] = binMakerSALT2_FORSALTREANALYSIS(All, neuron, stimBlock);
        [p I] = SALT(spt_baseline,spt_test,0.001,0.01);
        All{1,1}{neuron,1}(stimBlock,3) = p;
        All{1,1}{neuron,1}(stimBlock,4) = FSLatency;
        All{1,1}{neuron,1}(stimBlock,5) = jitter;
        All{1,1}{neuron,1}(stimBlock,6) = reliability;
        All{26,1}{neuron,2} = spt_spikeIx;
        BL_mean = mean(mean(spt_baseline,2),1);
        BL_std = std(mean(spt_baseline,2),1);
        T_peak = max(mean(spt_test(:,1:7),1));
        T_treshold = BL_mean+(1.65*BL_std);
        if T_peak > T_treshold
            T_Result = 1;
        elseif T_peak <= T_treshold
            T_Result = 0;
        end
        
        All{24,1}{neuron,1}(stimBlock,1) = BL_mean;
        All{24,1}{neuron,1}(stimBlock,2) = BL_std;
        All{24,1}{neuron,1}(stimBlock,3) = T_peak;
        All{24,1}{neuron,1}(stimBlock,4) = T_treshold;
        All{24,1}{neuron,1}(stimBlock,5) = T_Result;
        try
            [r1,p2] = corrcoef(All{13, 1}(neuron,(45:90)),All{16, 1}{neuron, stimBlock}(1,45:90));
             All{1,1}{neuron,1}(stimBlock,15) = r1(1,2);
        catch
            All{1,1}{neuron,1}(stimBlock,15) = nan;
        end
    end
end