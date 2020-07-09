function OIinfo = Lat_Jit_Rel_rangePlot(All, neuronIndex, parameters_OI)

% Function looks for OI neurons with sessions with 10 laser protocols with decreasing laser intensity. It chooses neurons
% that have more than one good OI protocol and has a very high chance of
% being identified based on very short lat, jit and Rel. Then plots the the latency, jitter and reliability for the maximal light vs. minimal light intensity where it was identified.
% Can be used to define your thresholds for optogenetic identification criteria

% Collect SALT info on sessions and laser protocols in sessions and store
% in matrix.
OIinfo = [];
m = 1;
for n = 1:size(All{1, 1},1)
if size(All{1, 1}{n, 1},1) == 10 && numel(find(ismember(n,neuronIndex.PyrIx_OI))) == 1
    lightIxs = find(All{1,1}{n,1}(:,3) < 0.01 & All{1,1}{n,1}(:,15)>= parameters_OI.PWCor);
    lightIxs = find(All{24,1}{n,1}(lightIxs,5) == 1);
    %lightIxs = lightIxs(1):lightIxs(end);
    if numel(lightIxs) > 3% & numel(find(All{1,1}{n,1}(:,4) < BestLatTreshold & All{1,1}{n,1}(:,5) < BestJitTreshold & All{1,1}{n,1}(:,6) > BestRelTreshold)) > 0
        [~, stimBlocksInt] = max(All{1,1}{n,1}(lightIxs,2));
        [~, stimBlocksDur] = max(All{1,1}{n,1}(lightIxs(stimBlocksInt),1));
        if m == 1
            OIinfo = All{1, 1}{n, 1}(lightIxs,1:6);
            OIinfo(stimBlocksInt(stimBlocksDur),7) = m;
            OIinfo(end,8) = m;
            OIinfo(:,9) = m;
            OIinfo(1,10) = n;
            m = m+1;
        else
            temp = All{1, 1}{n, 1}(lightIxs,1:6);
            temp(stimBlocksInt(stimBlocksDur),7) = m; % Indicate 'best' neuron with most light power and longest duration.
            temp(end,8) = m; % Indicate last protocol within sesions before light intesnity is too low.
            temp(:,9) = m; % Indicate # neuron it is in this matrix, can be used as unique identifier.
            temp(1,10) = n; % Indicate which neuron it was from in "All" struct.
            OIinfo(end+1:end+numel(lightIxs),:) = temp;
            m = m+1;
        end
        clear temp
    end
end
end

OIinfo([1 5 10 18 26 33 41 49],:) = [];
OIinfo(1,7) = 1;
OIinfo(4,7) = 2;
OIinfo(8,7) = 3;
OIinfo(15,7) = 4;
OIinfo(22,7) = 5;
OIinfo(28,7) = 6;
OIinfo(35,7) = 7;
OIinfo(42,7) = 8;


n=1;
for i=1:10
    Ix_last = find(OIinfo(:,8)==i);
    Ix_first = find(OIinfo(:,7)==i);
    Lat(i,1:2) = OIinfo([Ix_first,Ix_last],4);
    Jit(i,1:2) = OIinfo([Ix_first,Ix_last],5);
    Rel(i,1:2) = OIinfo([Ix_first,Ix_last],6);
    n = n+1;
end

figure
subplot(1,3,1)
hold on
plot(Lat(:,1:2)','k')
boxplot([Lat(:,1) Lat(:,2)],'Labels',{'Max light','Min light'})
ylabel('Latency (ms)')
ylim([0 7])
subplot(1,3,2)
hold on
plot(Jit(:,1:2)','k')
boxplot([Jit(:,1) Jit(:,2)],'Labels',{'Max light','Min light'})
ylabel('Jitter (ms)')
ylim([0 3.5])
subplot(1,3,3)
hold on
plot(Rel(:,1:2)','k')
boxplot([Rel(:,1) Rel(:,2)],'Labels',{'Max light','Min light'})
ylabel('Reliability')
ylim([0 1])

% Get Lat, Jit and Rel for 2 or 5 ms at 100 and 75 power
OIinfo = [];
m = 1;
for n = 1:size(All{1, 1},1)
if size(All{1, 1}{n, 1},1) >= 2 && numel(find(ismember(n,neuronIndex.PyrIx_OI))) == 1
    lightIxs = find(All{1,1}{n,1}(:,3) < 0.01 & All{1,1}{n,1}(:,15)>= parameters_OI.PWCor & All{1,1}{n,1}(:,1) == 2 & All{1,1}{n,1}(:,2) >= 75);
    %lightIxs = lightIxs(1):lightIxs(end);
    if numel(lightIxs) >= 2% & numel(find(All{1,1}{n,1}(:,4) < BestLatTreshold & All{1,1}{n,1}(:,5) < BestJitTreshold & All{1,1}{n,1}(:,6) > BestRelTreshold)) > 0
        if m == 1
            OIinfo = All{1, 1}{n, 1}(lightIxs,1:6);
            OIinfo(stimBlocksInt(stimBlocksDur),7) = m;
            OIinfo(end,8) = m;
            OIinfo(:,9) = m;
            OIinfo(1,10) = n;
            m = m+1;
        else
            temp = All{1, 1}{n, 1}(lightIxs,1:6);
            temp(stimBlocksInt(stimBlocksDur),7) = m; % Indicate 'best' neuron with most light power and longest duration.
            temp(end,8) = m; % Indicate last protocol within sesions before light intesnity is too low.
            temp(:,9) = m; % Indicate # neuron it is in this matrix, can be used as unique identifier.
            temp(1,10) = n; % Indicate which neuron it was from in "All" struct.
            OIinfo(end+1:end+numel(lightIxs),:) = temp;
            m = m+1;
        end
        clear temp
    end
end
end

n=1;
for i=1:OIinfo(end,9)
    Ix_last = find(OIinfo(:,8)==i);
    Ix_first = find(OIinfo(:,7)==i);
    Lat(i,1:2) = OIinfo([Ix_first,Ix_last],4);
    Jit(i,1:2) = OIinfo([Ix_first,Ix_last],5);
    Rel(i,1:2) = OIinfo([Ix_first,Ix_last],6);
    n = n+1;
end

figure
subplot(1,3,1)
hold on
for i=1:size(Lat,1)
    plot(Lat(i,1:2),'k','LineWidth',1)
end
plot(nanmedian(Lat(:,1:2),1),'r','LineWidth',2)
ylim([0 7])
ylabel('Latency(ms)')
subplot(1,3,2)
hold on
for i=1:size(Lat,1)
    plot(Jit(i,1:2),'k','LineWidth',1)
end
plot(nanmedian(Jit(:,1:2),1),'r','LineWidth',2)
ylim([0 3.5])
ylabel('Jitter(ms)')
subplot(1,3,3)
hold on
for i=1:size(Lat,1)
    plot(Rel(i,1:2),'k','LineWidth',1)
end
plot(nanmedian(Rel(:,1:2),1),'r','LineWidth',2)
ylim([0 1])
ylabel('Reliability')
