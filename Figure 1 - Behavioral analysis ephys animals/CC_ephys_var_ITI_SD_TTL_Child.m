function [Correct, Incorrect, Omission, Premature, Perserverative, nrTrialsStage] = CC_ephys_var_ITI_SD_TTL_Child(data,block)

% Function used by "CC_ephys_var_ITI_SD_TTL_Master". This functions
% processes each raw data array that the Master file provides and feeds
% back analysed data in structs per trial type.

pos = [2:5:length(data)];
ITI = unique(data(pos,1));
SD = unique(data(pos,5));
trialITI = data(pos,1);
trialSD = data(pos,5);

for l = 1:length(ITI)
    s = num2str(ITI(l)*10);
    ITIname = sprintf('TrialPosITI%s', s);
    Correct.(ITIname) = [];
    Incorrect.(ITIname) = [];
    Omission.(ITIname) = [];
    Premature.(ITIname) = [];
    Perserverative.(ITIname) = [];

    for i = pos
        varITIindex.(ITIname) = find(ismember(trialITI, ITI(l)));
    end
end

for l = 1:length(SD)
    s = num2str(SD(l)*10);
    SDname = sprintf('TrialPosSD%s', s);
    Correct.(SDname) = [];
    Incorrect.(SDname) = [];
    Omission.(SDname) = [];
    Premature.(SDname) = [];
    Perserverative.(SDname) = [];

    for i = pos
        varSDindex.(SDname) = find(ismember(trialSD, SD(l)));
    end
end

%% # trials for stage.
nrTrialsStage.TotalNr = length(pos);
%% Correct data for stage
pos = [3:5:length(data)];
Correct.TotalNr = length(find(ismember(data(pos,1), 2)));
Correct.TrialPosOverall = find(ismember(data(pos,1), 2));

n = 1;
for m = 1:length(Correct.TrialPosOverall)
    for l = 1:length(ITI)
    s = num2str(ITI(l)*10);
    ITIname = sprintf('TrialPosITI%s', s);
    %ITIname2 = sprintf('ITItimeIRbeam%s', s);
        if mean(ismember(varITIindex.(ITIname), Correct.TrialPosOverall(m)))>0
            Correct.(ITIname)(length(Correct.(ITIname))+1) = n;
            %Correct.(ITIname2)(length(Correct.(ITIname))+1) = data(pos(Correct.TrialPosOverall(m)),2);
            n = n+1;
        end
    end
end

n = 1;
for m = 1:length(Correct.TrialPosOverall)
    for l = 1:length(SD)
    s = num2str(SD(l)*10);
    SDname = sprintf('TrialPosSD%s', s);
    %SDname2 = sprintf('SDtimeIRbeam%s', s);
        if mean(ismember(varSDindex.(SDname), Correct.TrialPosOverall(m)))>0
            Correct.(SDname)(length(Correct.(SDname))+1) = n;
            %Correct.(SDname2)(length(Correct.(SDname))+1) = data(pos(Correct.TrialPosOverall(m)),2);
            n = n+1;
        end
    end
end
%% Incorrect data for stage
Incorrect.TotalNr = length(find(ismember(data(pos,1), 4)));
Incorrect.TrialPosOverall = find(ismember(data(pos,1), 4));

n = 1;
for m = 1:length(Incorrect.TrialPosOverall)
    for l = 1:length(ITI)
    s = num2str(ITI(l)*10);
    ITIname = sprintf('TrialPosITI%s', s);
    %ITIname2 = sprintf('ITItimeIRbeam%s', s);
        if mean(ismember(varITIindex.(ITIname), Incorrect.TrialPosOverall(m)))>0
            Incorrect.(ITIname)(length(Incorrect.(ITIname))+1) = n;
            %Incorrect.(ITIname2)(length(Incorrect.(ITIname))+1) = data(pos(Incorrect.TrialPosOverall(m)),2);
            n = n+1;
        end
    end
end

n = 1;
for m = 1:length(Incorrect.TrialPosOverall)
    for l = 1:length(SD)
    s = num2str(SD(l)*10);
    SDname = sprintf('TrialPosSD%s', s);
    %SDname2 = sprintf('SDtimeIRbeam%s', s);
        if mean(ismember(varSDindex.(SDname), Incorrect.TrialPosOverall(m)))>0
            Incorrect.(SDname)(length(Incorrect.(SDname))+1) = n;
            %Incorrect.(SDname2)(length(Incorrect.(SDname))+1) = data(pos(Incorrect.TrialPosOverall(m)),2);
            n = n+1;
        end
    end
end
%% Omission data for stage
Omission.TotalNr = length(find(ismember(data(pos,1), 1)));
Omission.TrialPosOverall = find(ismember(data(pos,1), 1));

n = 1;
for m = 1:length(Omission.TrialPosOverall)
    for l = 1:length(ITI)
    s = num2str(ITI(l)*10);
    ITIname = sprintf('TrialPosITI%s', s);
    %ITIname2 = sprintf('ITItimeIRbeam%s', s);
        if mean(ismember(varITIindex.(ITIname), Omission.TrialPosOverall(m)))>0
            Omission.(ITIname)(length(Omission.(ITIname))+1) = n;
            %Omission.(ITIname2)(length(Omission.(ITIname))+1) = data(pos(Omission.TrialPosOverall(m)),2);
            n = n+1;
        end
    end
end

n = 1;
for m = 1:length(Omission.TrialPosOverall)
    for l = 1:length(SD)
    s = num2str(SD(l)*10);
    SDname = sprintf('TrialPosSD%s', s);
    %SDname2 = sprintf('SDtimeIRbeam%s', s);
        if mean(ismember(varSDindex.(SDname), Omission.TrialPosOverall(m)))>0
            Omission.(SDname)(length(Omission.(SDname))+1) = n;
            %Omission.(SDname2)(length(Omission.(SDname))+1) = data(pos(Omission.TrialPosOverall(m)),2);
            n = n+1;
        end
    end
end
%% Premature response data for stage
Premature.TotalNr = length(find(ismember(data(pos,1), 3)));
Premature.TrialPosOverall = find(ismember(data(pos,1), 3));

n = 1;
for m = 1:length(Premature.TrialPosOverall)
    for l = 1:length(ITI)
    s = num2str(ITI(l)*10);
    ITIname = sprintf('TrialPosITI%s', s);
    %ITIname2 = sprintf('ITItimeIRbeam%s', s);
        if mean(ismember(varITIindex.(ITIname), Premature.TrialPosOverall(m)))>0
            Premature.(ITIname)(length(Premature.(ITIname))+1) = n;
            %Premature.(ITIname2)(length(Premature.(ITIname))+1) = data(pos(Premature.TrialPosOverall(m)),2);
            n = n+1;
        end
    end
end

n = 1;
for m = 1:length(Premature.TrialPosOverall)
    for l = 1:length(SD)
    s = num2str(SD(l)*10);
    SDname = sprintf('TrialPosSD%s', s);
    %SDname2 = sprintf('SDtimeIRbeam%s', s);
        if mean(ismember(varSDindex.(SDname), Premature.TrialPosOverall(m)))>0
            Premature.(SDname)(length(Premature.(SDname))+1) = n;
            %Premature.(SDname2)(length(Premature.(SDname))+1) = data(pos(Premature.TrialPosOverall(m)),2);
            n = n+1;
        end
    end
end
%% Perserverative data for stage
pos = [5:5:length(data)];
Perserverative.TrialPosOverall = find(data(pos,:));
n = 1;
for i = 1:length(pos)
    if sum(data(pos(i),:))>0
        Perserverative.TrialPosOverall(n) = i;
        n = n+1;
    end
end
Perserverative.TotalNr = length(Perserverative.TrialPosOverall);

n = 1;
for m = 1:length(Perserverative.TrialPosOverall)
    for l = 1:length(ITI)
    s = num2str(ITI(l)*10);
    ITIname = sprintf('TrialPosITI%s', s);
    %ITIname2 = sprintf('ITItimeIRbeam%s', s);
        if mean(ismember(varITIindex.(ITIname), Perserverative.TrialPosOverall(m)))>0
            Perserverative.(ITIname)(length(Perserverative.(ITIname))+1) = n;
            %Perserverative.(ITIname2)(length(Perserverative.(ITIname))+1) = data(pos(Perserverative.TrialPosOverall(m)),2);
            n = n+1;
        end
    end
end

n = 1;
for m = 1:length(Perserverative.TrialPosOverall)
    for l = 1:length(SD)
    s = num2str(SD(l)*10);
    SDname = sprintf('TrialPosSD%s', s);
    %SDname2 = sprintf('SDtimeIRbeam%s', s);
        if mean(ismember(varSDindex.(SDname), Perserverative.TrialPosOverall(m)))>0
            Perserverative.(SDname)(length(Perserverative.(SDname))+1) = n;
            %Perserverative.(SDname2)(length(Perserverative.(SDname))+1) = data(pos(Perserverative.TrialPosOverall(m)),2);
            n = n+1;
        end
    end
end

%% Get correct latency times
n = 1;
for i = 1:numel(Correct.TrialPosOverall)
    pos = [2:5:length(data)];
    Correct.MagLatency(n) = data(pos(Correct.TrialPosOverall(n)),4);
    n = n+1;
end

conditions = fieldnames(Correct);
check = 1;

for i = 1:numel(conditions)
    if strcmp(conditions{i}(end-5:end),'ITI125')
        check = 1;
        break
    elseif strcmp(conditions{i}(end-2:end),'SD2')
        check = 0;
        break
    else
        check = 2;
    end
end

if check == 1
    Correct.MagLatencyITI50 = Correct.MagLatency(Correct.TrialPosITI50);
    Correct.MagLatencyITI75 = Correct.MagLatency(Correct.TrialPosITI75);
    Correct.MagLatencyITI125 = Correct.MagLatency(Correct.TrialPosITI125);
    Correct.MagLatencySD10 = Correct.MagLatency(Correct.TrialPosSD10);
elseif check == 0
    Correct.MagLatencySD10 = Correct.MagLatency(Correct.TrialPosSD10);
    Correct.MagLatencySD5 = Correct.MagLatency(Correct.TrialPosSD5);
    Correct.MagLatencySD2 = Correct.MagLatency(Correct.TrialPosSD2);
    Correct.MagLatencyITI5 = Correct.MagLatency(Correct.TrialPosITI50);
end

%% Get premature latency times
n = 1;
for i = 1:numel(Premature.TrialPosOverall)
    pos_ETO = [4:5:length(data)];
    pos_STr = [1:5:length(data)];
    Premature.Latency(n) = data(pos_ETO(Premature.TrialPosOverall(n)),2)-data(pos_STr(Premature.TrialPosOverall(n)),1)-5;
    n = n+1;
end

conditions = fieldnames(Premature);
check = 1;

for i = 1:numel(conditions)
    if strcmp(conditions{i}(end-5:end),'ITI125')
        check = 1;
        break
    elseif strcmp(conditions{i}(end-2:end),'SD2')
        check = 0;
        break
    else
        check = 2;
    end
end

if check == 1
    Premature.LatencyITI50 = Premature.Latency(Premature.TrialPosITI50);
    Premature.LatencyITI75 = Premature.Latency(Premature.TrialPosITI75);
    Premature.LatencyITI125 = Premature.Latency(Premature.TrialPosITI125);
    Premature.LatencySD10 = Premature.Latency(Premature.TrialPosSD10);
elseif check == 0
    Premature.LatencySD10 = Premature.Latency(Premature.TrialPosSD10);
    Premature.LatencySD5 = Premature.Latency(Premature.TrialPosSD5);
    Premature.LatencySD2 = Premature.Latency(Premature.TrialPosSD2);
    Premature.LatencyITI5 = Premature.Latency(Premature.TrialPosITI50);
end
