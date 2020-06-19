function [spt_test, spt_baseline, FSLatency, jitter, reliability, spt_spikeIx, baselineTimeArray] = binMakerSALT2(All, eventArray, spikes)


%RecParam;

% if stimBlock == 1;
%     range = stimBlock:500;
% else
%     range = 1+((stimBlock-1)*500):((stimBlock-1)*500)+500;
% end

% Get event timestamps
% timeStamps = get_eventTimesNeurorighter;
% CH = 'CH16';
% eventArray = All{22, 1}{clusterNumber,1}.(CH)(range);
baselineTimeArrayTEMP = eventArray;
baselineTimeArray = eventArray;

randBaselineOffset = linspace(0.02, 0.1, 20);

for i = 1:length(baselineTimeArrayTEMP)
    randTime = datasample(randBaselineOffset,1);
    baselineTimeArray(i) = baselineTimeArrayTEMP(i)-0.04;
end

% Get spike timestamps
% spikes = All{21, 1}{clusterNumber,1};

% Define parameters
binSize = 0.001;
segmentSize = 20;

% Pre-allocate test and baseline matrix that are input for SALT test.
spt_test = zeros(length(eventArray),segmentSize);
spt_baseline = zeros(length(eventArray),segmentSize*2);

spt_spikeIx = [];
n=1;
p=1;
FSLatencyAll = zeros(1,length(eventArray));
% Binarize spikes after event timestamps and store in spt_test
for i = 1:length(eventArray)
    % Skip light off artefact bin, depending on stim length
     j = linspace(1,segmentSize,segmentSize);
%     if Stim.Dur(stimBlock) < segmentSize
%         j(Stim.Dur(stimBlock)+1) = [];
%     end
    for m = 1:length(j) %skip first millisecond for light on artefact
        x = find(spikes>=(eventArray(i)+((j(m)-1)*binSize)) & spikes<eventArray(i)+((j(m)-1)*binSize)+binSize);
        if ~isempty(x)
            spt_test(i,j(m))=1;
            n = n+1;
            if j(m)<=10
                % Store first spike latency
                FSLatencyAll(1,i) = j(m);
                spt_spikeIx(1,p) = x(1);
                p = p+1;
            end
            break
        end
    end
end

% Binarize spikes after event timestamps and store in spt_baseline
for i = 1:length(baselineTimeArray)
    for l = 1:(segmentSize*2)
        x = find(spikes>=(baselineTimeArray(i)+((l-1)*binSize)) & spikes<baselineTimeArray(i)+((l-1)*binSize)+binSize);
        if ~isempty(x)
            spt_baseline(i,l)=1;
            break
        end
    end
end

% Arteficially increase trials by random resampling and arteficially
% increase length of baseline array by random resampling
% i = 1;
% while i < 100
%     spt_baseline = [spt_baseline; spt_baseline(datasample(1:500,500),:)];
%     i = i+1;
% end

spt_baseline = [spt_baseline zeros(length(spt_baseline(:,1)),760)];

i = 1;
while i < 20
    spt_baseline(:,i*40+1:(i+1)*40) = spt_baseline(:,datasample(1:40,40));
    i = i+1;
end
    

% Latency of first spikes in ms.
FSLatency = sum(FSLatencyAll)/length(find(FSLatencyAll>0));
% Jitter calculated as standard deviations in spike times (zeros exluded).
jitter = std(FSLatencyAll(find(FSLatencyAll>0)));
% First spike reliability #first spikes/total light stims.
reliability = length(find(FSLatencyAll>0))/length(FSLatencyAll);