function [neuronIndex, parameters_OI, All] = Cluster_sort(All)
% Apply cluster quality metrics per neuron, optogenetic identification criteria per neuron and seperate in broad-spiking, narrow-spiking and unclassified.
% ALso puts distance from pia per neuron in master struct "All"

%% Apply cluster quality filter based on L-Ratio treshold and ISI interval violation
neuronIndex.ClustIx = 1:size(All{1,1},1);
neuronIndex.ClustIxBAD = [261 389 405 427 835 861 874 876 881 895 904 907 914 949 959 963 964 974 1003 1045 1060 1088 1103 1205 1217 1221 1225 1252 1254 1255 1260 1262 1264 1267 1278 1282 1372 1385 1438 1452 1453 1467,...
    4 71 72 76 88 135 242 254 446 473 480 492 498 533 558 563 565 582 593 657 667 679 691 753 810 1120 1123 1125 1144 1146 1147 1174 1176 1200 1298 1316 1343 1345];

% Indicate OI criteria
parameters_OI.PWCor = 0.80; %Minimal pearon's waveform correlation
parameters_OI.minLightInt = 75; %Minimal laser power requried for optogenetic identification
parameters_OI.minSessionSALT = 1; %Minimal significant SALT sessions per neuron before being accepted
parameters_OI.JitterTreshold = 3.5; %Maximal first spike jitter threshold
parameters_OI.LatencyTreshold = 7; %Maximal first spike altency threshold
parameters_OI.ReliabilityTreshold = 0.031; %low firing units in baseline period tend to give false positives. Manual inspection seem to indicate...
%that Reliability of 0.031 (~15 spikes in 500 trials) is a good lower
%treshold.

neuronIndex.manualOIinspection = []; %manually filter out false positive OI neurons

neuronIndex.manualOIinspectionInsert = [];

% Indicate cluster criteria
neuronIndex.LRtreshold = 1.5;
neuronIndex.IDtreshold = 40;
neuronIndex.ISI15treshold = 0.25; %Percentage of spikes in ISI <1.5ms

%% END OF PARAMETERS

for i = 1:size(All{1,1},1)
    if All{3,1}(i,2) > neuronIndex.LRtreshold || All{3,1}(i,7) > neuronIndex.ISI15treshold || All{3,1}(i,1) < neuronIndex.IDtreshold
        neuronIndex.ClustIxBAD = [neuronIndex.ClustIxBAD i];
    end
end

neuronIndex.ClustIx(neuronIndex.ClustIxBAD) = [];
neuronIndex.Clust_All_reliable = neuronIndex.ClustIx;

        
%% Find OI neurons based on the following: minimally 75% light power, SALT p<0.01 and pearson's waveform corr. >=0.80
neuronIndex.OIindex = [];
neuronIndex.OIindexMaybe = [];

for i = neuronIndex.ClustIx 
    lightIxs = find(All{1,1}{i,1}(:,3) < 0.01 & All{1,1}{i,1}(:,15)>= parameters_OI.PWCor & All{1,1}{i,1}(:,2)>= parameters_OI.minLightInt & All{1,1}{i,1}(:,6)>= parameters_OI.ReliabilityTreshold & All{1,1}{i,1}(:,4) < parameters_OI.LatencyTreshold & All{1,1}{i,1}(:,5) < parameters_OI.JitterTreshold);
    lightIxsMaybe = find(All{1,1}{i,1}(:,3) < 0.05 & All{1,1}{i,1}(:,3) >= 0.01 & All{1,1}{i,1}(:,15)>= parameters_OI.PWCor & All{1,1}{i,1}(:,2)>= parameters_OI.minLightInt & All{1,1}{i,1}(:,6)>= parameters_OI.ReliabilityTreshold & All{1,1}{i,1}(:,4) < parameters_OI.LatencyTreshold & All{1,1}{i,1}(:,5) < parameters_OI.JitterTreshold);

    if numel(lightIxs) >= parameters_OI.minSessionSALT && numel(find(All{24,1}{i,1}(lightIxs,5) == 1)) > 0
        neuronIndex.OIindex = [neuronIndex.OIindex i];
    end
    if numel(lightIxsMaybe) >= parameters_OI.minSessionSALT
        neuronIndex.OIindexMaybe = [neuronIndex.OIindexMaybe i];
    end
end

for i = 1:length(neuronIndex.OIindex)
    neuronIndex.ClustIx(find(neuronIndex.ClustIx==neuronIndex.OIindex(i))) = [];
end

neuronIndex.ClustIx = [neuronIndex.ClustIx neuronIndex.manualOIinspection];
 
for i = neuronIndex.manualOIinspection
    neuronIndex.OIindex(find(neuronIndex.OIindex==i)) = [];
end

%% Find cluster indexes for each behavioral condition (vITI, vSD and FIXED_ITI) for cluster with and without OI
neuronIndex.vITIindex = [];
neuronIndex.vSDindex = [];

neuronIndex.vITIindexOI = [];
neuronIndex.vSDindexOI = [];
for i = neuronIndex.ClustIx
    if strcmp(All{5,1}{i,1},'VAR_ITI');
        neuronIndex.vITIindex = [neuronIndex.vITIindex i];
    elseif strcmp(All{5,1}{i,1},'VAR_SD');
        neuronIndex.vSDindex = [neuronIndex.vSDindex i];
    end
end

for i = neuronIndex.OIindex
    if strcmp(All{5,1}{i,1},'VAR_ITI');
        neuronIndex.vITIindexOI = [neuronIndex.vITIindexOI i];
    elseif strcmp(All{5,1}{i,1},'VAR_SD');
        neuronIndex.vSDindexOI = [neuronIndex.vSDindexOI i];
    end
end

%% Recalculate mean firing rate based on only behavior
for neuron = 1:size(All{15, 1},1)
All{15, 1}(neuron,2) = All{3, 1}(neuron,3)/All{21, 1}{neuron, 1}(All{3, 1}(neuron,3));
end

%% Assign putative pyramidal or interneuron with GMM (all neurons all sessions) clustering ans store in 23th cell as '1' - pyramidal or '2' - interneuron
X = [All{15,1}(neuronIndex.Clust_All_reliable,5),All{15,1}(neuronIndex.Clust_All_reliable,6)];

% Fit Gaussian mixture model and look for 2 clusters
k =2;
%Sigma = {'diagonal', 'full'};
Sigma = {'diagonal'};
nSigma = numel(Sigma);
% SharedCovariance = {true, false};
% SCtext = {'true', 'false'};
SharedCovariance = {false};
SCtext = {'false'};
nSC = numel(SharedCovariance);
d = 500;
x1 = linspace(min(X(:,1)) - 2,max(X(:,1)) + 2,d);
x2 = linspace(min(X(:,2)) - 2,max(X(:,2)) + 2,d);
[x1grid,x2grid] = meshgrid(x1,x2);
X0 = [x1grid(:) x2grid(:)];
threshold = sqrt(chi2inv(0.99,2));
options = statset('MaxIter',1000); % Increase number of EM iterations
%figure;
%c = 1;
for i = 1:nSigma;
    for j = 1:nSC;
        gmfit = fitgmdist(X,k,'CovarianceType',Sigma{i},...
        'SharedCovariance',SharedCovariance{j},'Options',options,'Replicates',10);
        clusterX = cluster(gmfit,X);
        mahalDist = mahal(gmfit,X0);
        for m = 1:k;
            idx = mahalDist(:,m)<=threshold;
        end
    end
end
P = posterior(gmfit,X);

% Remove any units with classification confidence <0.85 and classify as
% type 3 neuron
for i = 1:size(P,1)
    if max(P(i,:))<0.85
        clusterX(i) = 3;
    end
end

% Look for type 1 and 2 neurons and put in "All" struct
All{23,1} = nan(size(All{22, 1},1),1);
All{23,2} = 'Pyramidal (1), interneuron (2) or undefined (3)';

NrOnes = size(find(clusterX==1),1);
NrTwos = size(find(clusterX==2),1);
[z, pyrNr] = max([NrOnes NrTwos]);
[z, intNr] = min([NrOnes NrTwos]);


neuronPos = neuronIndex.Clust_All_reliable';
putneuronIndex.PyrIx = find(clusterX==pyrNr);
putneuronIndex.IntIx = find(clusterX==intNr);
putUnIDIx = find(clusterX==3);

%NrOnes = size(find(idx==1),1);
%NrTwos = size(find(idx==2),1);
%[z, pyrNr] = max([NrOnes NrTwos]);
%[z, intNr] = min([NrOnes NrTwos]);
All{23,1}(neuronPos(putneuronIndex.PyrIx)) = 1;
All{23,1}(neuronPos(putneuronIndex.IntIx)) = 2;
All{23,1}(neuronPos(putUnIDIx),1) = 3;

% Sort neurons on putative pyramidal or interneuron

neuronIndex.ClustIxPyr = find(All{23,1} == 1);
neuronIndex.ClustIxInt = find(All{23,1} == 2);
neuronIndex.ClustIxUnknown = find(All{23,1} == 3);

neuronIndex.PyrIx_OI = neuronIndex.ClustIxPyr(find(ismember(neuronIndex.ClustIxPyr, neuronIndex.OIindex)));
neuronIndex.IntIx_OI = neuronIndex.ClustIxInt(find(ismember(neuronIndex.ClustIxInt, neuronIndex.OIindex)));
neuronIndex.PyrIx = neuronIndex.ClustIxPyr(find(ismember(neuronIndex.ClustIxPyr, neuronIndex.ClustIx)));
neuronIndex.IntIx = neuronIndex.ClustIxInt(find(ismember(neuronIndex.ClustIxInt, neuronIndex.ClustIx)));

neuronIndex.PyrIx_VAR_ITI_OI = neuronIndex.vITIindexOI(find(ismember(neuronIndex.vITIindexOI, neuronIndex.PyrIx_OI)));
neuronIndex.PyrIx_VAR_SD_OI = neuronIndex.vSDindexOI(find(ismember(neuronIndex.vSDindexOI, neuronIndex.PyrIx_OI)));
neuronIndex.IntIxVAR_ITI_OI = neuronIndex.vITIindexOI(find(ismember(neuronIndex.vITIindexOI, neuronIndex.IntIx_OI)));
neuronIndex.IntIxVAR_SD_OI = neuronIndex.vSDindexOI(find(ismember(neuronIndex.vSDindexOI, neuronIndex.IntIx_OI)));
neuronIndex.PyrIxVAR_ITI = neuronIndex.vITIindex(find(ismember(neuronIndex.vITIindex, neuronIndex.PyrIx)));
neuronIndex.PyrIxVAR_SD = neuronIndex.vSDindex(find(ismember(neuronIndex.vSDindex, neuronIndex.PyrIx)));
neuronIndex.IntIxVAR_ITI = neuronIndex.vITIindex(find(ismember(neuronIndex.vITIindex, neuronIndex.IntIx)));
neuronIndex.IntIxVAR_SD = neuronIndex.vSDindex(find(ismember(neuronIndex.vSDindex, neuronIndex.IntIx)));

%lat_jit_Rel_KTEST2 = waveformAnalysis(All, parameters_OI, neuronIndex, neuronIndex.IntIxVAR_ITI, neuronIndex.IntIxVAR_SD, neuronIndex.PyrIxVAR_ITI, neuronIndex.PyrIxVAR_SD, neuronIndex.PyrIx_VAR_ITI_OI, neuronIndex.PyrIx_VAR_SD_OI, neuronIndex.IntIxVAR_ITI_OI, neuronIndex.IntIxVAR_SD_OI);

% Plot plot cluster quality of all(none excluded yet) neurons and OI
% neurons to show that cluster criteria are strict
neuronIndex.OIindexTEMP = [];

for i = 1:size(All{1,1},1)
    lightIxs = find(All{1,1}{i,1}(:,3) < 0.01 & All{1,1}{i,1}(:,15)>= parameters_OI.PWCor & All{1,1}{i,1}(:,2)>= parameters_OI.minLightInt & All{1,1}{i,1}(:,6)>= parameters_OI.ReliabilityTreshold);
    if numel(lightIxs) >= parameters_OI.minSessionSALT
        neuronIndex.OIindexTEMP = [neuronIndex.OIindexTEMP i];
    end
end

%% Find pairs of neurons that are possibly the same across VAR_ITI and VAR_SD sessions
for neuron = 1:size(neuronIndex.ClustIxPyr,1)
    AnimalNeurons = All{8, 1}{neuronIndex.ClustIxPyr(neuron),1};
    DepthNeurons = All{11, 1}(neuronIndex.ClustIxPyr(neuron),1);
    shankNeurons = All{4, 1}(neuronIndex.ClustIxPyr(neuron),1);
    dateNeurons = All{9, 1}{neuronIndex.ClustIxPyr(neuron),1};
    protNeurons = All{5, 1}{neuronIndex.ClustIxPyr(neuron), 1};
    peakCHNeurons = All{15, 1}(neuronIndex.ClustIxPyr(neuron), 1);
    
    targetNeurons = [];
    for i = 1:size(neuronIndex.ClustIxPyr,1)
        if All{8, 1}{neuronIndex.ClustIxPyr(i),1} == AnimalNeurons & All{11, 1}(neuronIndex.ClustIxPyr(i),1) == DepthNeurons & All{4, 1}(neuronIndex.ClustIxPyr(i),1) == shankNeurons & strcmp(All{9, 1}{neuronIndex.ClustIxPyr(i),1},dateNeurons) == 0 & strcmp(All{5, 1}{neuronIndex.ClustIxPyr(i),1},protNeurons) == 0 & (All{15,1}(neuronIndex.ClustIxPyr(i),1)>=peakCHNeurons-1 & All{15,1}(neuronIndex.ClustIxPyr(i),1)<=peakCHNeurons+1)
            targetNeurons = [targetNeurons neuronIndex.ClustIxPyr(i)];
        end
    end
    
    if ~isempty(targetNeurons)
        WVcor = [];
        for partnerNeuron = 1:numel(targetNeurons)
            R = corrcoef(All{13, 1}(neuronIndex.ClustIxPyr(neuron),:),All{13, 1}(targetNeurons(partnerNeuron),:));
            WVcor(partnerNeuron) = R(1,2);
        end

        WvCorIx = find(WVcor>0.95);
        if ~isempty(WvCorIx)
            ISIcor = [];
            for i = 1:numel(WvCorIx)
                spikeIntervalsTarget = diff(All{21,1}{targetNeurons(WvCorIx(i)),1}(1:All{3, 1}(targetNeurons(WvCorIx(i)),3)))*1000;
                spikeIntervalsNeuron = diff(All{21,1}{neuronIndex.ClustIxPyr(neuron),1}(1:All{3, 1}(neuronIndex.ClustIxPyr(neuron),3)))*1000;
                binSize_hist = 0.5;
                x = [0:binSize_hist:100];
                targetNeuronHist = histcounts(spikeIntervalsTarget,x);
                NeuronHist = histcounts(spikeIntervalsNeuron,x);
                R = corrcoef(targetNeuronHist, NeuronHist);
                ISIcor(i) = R(1,2);
            end

            if max(ISIcor)>0.70
                [C I] = max(ISIcor);
                %All{25,1}(neuronIndex.ClustIxPyr(neuron),:) = targetNeurons(WvCorIx(I));
            else
                All{25,1}(neuronIndex.ClustIxPyr(neuron),:) = nan;
            end

            if contains(protNeurons,'VAR_SD') && max(ISIcor)>0.70
                [~, PSTH_neuron, ~] = PSTH_event_SD(All, neuronIndex, parameters_OI, neuronIndex.ClustIxPyr(neuron), 0.2, 'All', 'no');
                [~, PSTH_targetNeuron, ~] = PSTH_event_ITI(All, neuronIndex, parameters_OI, targetNeurons(WvCorIx(I)), 0.2, 'All', 'no');
                R_neuron = corrcoef(PSTH_neuron.Cor.wait_start(6:25), PSTH_targetNeuron.Cor.wait_start(6:25));
                if R_neuron(1,2)>0.5
                    All{25,1}(neuronIndex.ClustIxPyr(neuron),:) = targetNeurons(WvCorIx(I));
                else
                    All{25,1}(neuronIndex.ClustIxPyr(neuron),:) = nan;
                end

            elseif contains(protNeurons,'VAR_ITI') && max(ISIcor)>0.70
                [~, PSTH_neuron, ~] = PSTH_event_ITI(All, neuronIndex, parameters_OI, neuronIndex.ClustIxPyr(neuron), 0.2, 'All', 'no');
                [~, PSTH_targetNeuron, ~] = PSTH_event_SD(All, neuronIndex, parameters_OI, targetNeurons(WvCorIx(I)), 0.2, 'All', 'no');
                R_neuron = corrcoef(PSTH_neuron.Cor.wait_start(6:25), PSTH_targetNeuron.Cor.wait_start(6:25));
                if R_neuron(1,2)>0.5
                    All{25,1}(neuronIndex.ClustIxPyr(neuron),:) = targetNeurons(WvCorIx(I));
                else
                    All{25,1}(neuronIndex.ClustIxPyr(neuron),:) = nan;
                end

            end
        end
    elseif isempty(partnerNeuron)
        All{25,1}(neuronIndex.ClustIxPyr(neuron),:) = nan;
    end
    clear targetNeurons parterNeuron R R_neuron
end

All{25,1}(All{25,1} == 0) = nan;

% Search for partner neurons that were selected twice. Run correlation over
% PSTH around treshold in that case.
% n = [];
% for neuron = 1:size(neuronIndex.ClustIxPyr,1)
%     if ~isnan(All{25,1}(neuronIndex.ClustIxPyr(neuron))) && numel(find(All{25,1} == All{25,1}(neuronIndex.ClustIxPyr(neuron))))>1
%         DoubleIxs = find(All{25,1} == All{25,1}(neuronIndex.ClustIxPyr(neuron)));
%         targetNeuron = All{25,1}(DoubleIxs(1));
%         protNeurons = All{5, 1}{neuronIndex.ClustIxPyr(neuron), 1};
%         if contains(protNeurons,'VAR_SD')
%             [~, PSTH_neuron, ~] = PSTH_event_SD(All, neuronIndex, parameters_OI, DoubleIxs(1), 0.2, 'All', 'no');
%             [~, PSTH_neuron_two, ~] = PSTH_event_SD(All, neuronIndex, parameters_OI, DoubleIxs(2), 0.2, 'All', 'no');
%             [~, PSTH_targetNeuron, ~] = PSTH_event_ITI(All, neuronIndex, parameters_OI, targetNeuron, 0.2, 'All', 'no');
%             R_neuron = corrcoef(PSTH_neuron.Cor.wait_start(6:25), PSTH_targetNeuron.Cor.wait_start(6:25));
%             R_neuron_two = corrcoef(PSTH_neuron_two.Cor.wait_start(6:25), PSTH_targetNeuron.Cor.wait_start(6:25));
%             [~, I] = max([R_neuron(1,2) R_neuron_two(1,2)]);
%             DoubleIxs(I) = [];
%             All{25,1}(DoubleIxs,1) = nan;
%         elseif contains(protNeurons,'VAR_ITI')
%             [~, PSTH_neuron, ~] = PSTH_event_ITI(All, neuronIndex, parameters_OI, DoubleIxs(1), 0.2, 'All', 'no');
%             [~, PSTH_neuron_two, ~] = PSTH_event_ITI(All, neuronIndex, parameters_OI, DoubleIxs(2), 0.2, 'All', 'no');
%             [~, PSTH_targetNeuron, ~] = PSTH_event_SD(All, neuronIndex, parameters_OI, targetNeuron, 0.2, 'All', 'no');
%             R_neuron = corrcoef(PSTH_neuron.Cor.wait_start(6:25), PSTH_targetNeuron.Cor.wait_start(6:25));
%             R_neuron_two = corrcoef(PSTH_neuron_two.Cor.wait_start(6:25), PSTH_targetNeuron.Cor.wait_start(6:25));
%             [~, I] = max([R_neuron(1,2) R_neuron_two(1,2)]);
%             DoubleIxs(I) = [];
%             All{25,1}(DoubleIxs,1) = nan;
%         end
%     end
% end

%% Insert distance from pia measure per neuron

for i = 1:size(All{4, 1},1)
    if contains(All{8, 1}{i,1},'DS-OI-2') & All{4, 1}(i,1) == 0
        All{4, 1}(i,2) = 616;
    elseif contains(All{8, 1}{i,1},'DS-OI-2') & All{4, 1}(i,1) == 1
        All{4, 1}(i,2) = 616+250;
    elseif contains(All{8, 1}{i,1},'DS-OI-2') & All{4, 1}(i,1) == 2
        All{4, 1}(i,2) = 616+500;
    elseif contains(All{8, 1}{i,1},'DS-OI-2') & All{4, 1}(i,1) == 3
        All{4, 1}(i,2) = 616+750;
        
    elseif contains(All{8, 1}{i,1},'DS-OI-5') & All{4, 1}(i,1) == 0
        All{4, 1}(i,2) = 476;
    elseif contains(All{8, 1}{i,1},'DS-OI-5') & All{4, 1}(i,1) == 1
        All{4, 1}(i,2) = 476+250;
    elseif contains(All{8, 1}{i,1},'DS-OI-5') & All{4, 1}(i,1) == 2
        All{4, 1}(i,2) = 476+500;
    elseif contains(All{8, 1}{i,1},'DS-OI-5') & All{4, 1}(i,1) == 3
        All{4, 1}(i,2) = 476+750;
        
    elseif contains(All{8, 1}{i,1},'DS-OI-6') & All{4, 1}(i,1) == 0
        All{4, 1}(i,2) = 362+750;
    elseif contains(All{8, 1}{i,1},'DS-OI-6') & All{4, 1}(i,1) == 1
        All{4, 1}(i,2) = 362+500;
    elseif contains(All{8, 1}{i,1},'DS-OI-6') & All{4, 1}(i,1) == 2
        All{4, 1}(i,2) = 362+250;
    elseif contains(All{8, 1}{i,1},'DS-OI-6') & All{4, 1}(i,1) == 3
        All{4, 1}(i,2) = 362;
        
    elseif contains(All{8, 1}{i,1},'DS-OI-7') & All{4, 1}(i,1) == 0
        All{4, 1}(i,2) = 612;
    elseif contains(All{8, 1}{i,1},'DS-OI-7') & All{4, 1}(i,1) == 1
        All{4, 1}(i,2) = 612+250;
    elseif contains(All{8, 1}{i,1},'DS-OI-7') & All{4, 1}(i,1) == 2
        All{4, 1}(i,2) = 612+500;
    elseif contains(All{8, 1}{i,1},'DS-OI-7') & All{4, 1}(i,1) == 3
        All{4, 1}(i,2) = 612+750;
        
    elseif contains(All{8, 1}{i,1},'DS-OI-8') & All{4, 1}(i,1) == 0
        All{4, 1}(i,2) = 802+750;
    elseif contains(All{8, 1}{i,1},'DS-OI-8') & All{4, 1}(i,1) == 1
        All{4, 1}(i,2) = 802+500;
    elseif contains(All{8, 1}{i,1},'DS-OI-8') & All{4, 1}(i,1) == 2
        All{4, 1}(i,2) = 802+250;
    elseif contains(All{8, 1}{i,1},'DS-OI-8') & All{4, 1}(i,1) == 3
        All{4, 1}(i,2) = 802;
    end
end

% Go over partnerneurons from OI neurons. If they are not marked as OI, put
% them there irrespective of SALT result. We assume here that OI is not
% perfect over days.

% for neuron = 1:numel(neuronIndex.PyrIx_VAR_SD_OI)
%     if ~isnan(All{25, 1}(neuronIndex.PyrIx_VAR_SD_OI(neuron),1)) && ~ismember(All{25, 1}(neuronIndex.PyrIx_VAR_SD_OI(neuron),1),neuronIndex.PyrIx_VAR_ITI_OI)
%         neuronIndex.PyrIx_VAR_ITI_OI = [neuronIndex.PyrIx_VAR_ITI_OI, All{25, 1}(neuronIndex.PyrIx_VAR_SD_OI(neuron),1)];
%         neuronIndex.PyrIx_OI = [neuronIndex.PyrIx_OI; All{25, 1}(neuronIndex.PyrIx_VAR_SD_OI(neuron),1)];
%         neuronIndex.OIindex = [neuronIndex.OIindex, All{25, 1}(neuronIndex.PyrIx_VAR_SD_OI(neuron),1)];
%         neuronIndex.vITIindexOI = [neuronIndex.vITIindexOI, All{25, 1}(neuronIndex.PyrIx_VAR_SD_OI(neuron),1)];
%         neuronIndex.PyrIx(find(neuronIndex.PyrIx==All{25, 1}(neuronIndex.PyrIx_VAR_SD_OI(neuron),1))) = [];
%         neuronIndex.PyrIxVAR_ITI(find(neuronIndex.PyrIxVAR_ITI==All{25, 1}(neuronIndex.PyrIx_VAR_SD_OI(neuron),1))) = [];
%         neuronIndex.vITIindex(find(neuronIndex.vITIindex==All{25, 1}(neuronIndex.PyrIx_VAR_SD_OI(neuron),1))) = [];
%     end
% end
% 
% for neuron = 1:numel(neuronIndex.PyrIx_VAR_ITI_OI)
%     if ~isnan(All{25, 1}(neuronIndex.PyrIx_VAR_ITI_OI(neuron),1)) && ~ismember(All{25, 1}(neuronIndex.PyrIx_VAR_ITI_OI(neuron),1),neuronIndex.PyrIx_VAR_SD_OI)
%         neuronIndex.PyrIx_VAR_SD_OI = [neuronIndex.PyrIx_VAR_SD_OI, All{25, 1}(neuronIndex.PyrIx_VAR_ITI_OI(neuron),1)];
%         neuronIndex.PyrIx_OI = [neuronIndex.PyrIx_OI; All{25, 1}(neuronIndex.PyrIx_VAR_ITI_OI(neuron),1)];
%         neuronIndex.OIindex = [neuronIndex.OIindex, All{25, 1}(neuronIndex.PyrIx_VAR_ITI_OI(neuron),1)];
%         neuronIndex.vSDindexOI = [neuronIndex.vSDindexOI, All{25, 1}(neuronIndex.PyrIx_VAR_ITI_OI(neuron),1)];
%         neuronIndex.PyrIx(find(neuronIndex.PyrIx==All{25, 1}(neuronIndex.PyrIx_VAR_ITI_OI(neuron),1))) = [];
%         neuronIndex.PyrIxVAR_SD(find(neuronIndex.PyrIxVAR_SD==All{25, 1}(neuronIndex.PyrIx_VAR_ITI_OI(neuron),1))) = [];
%         neuronIndex.vSDindex(find(neuronIndex.vSDindex==All{25, 1}(neuronIndex.PyrIx_VAR_ITI_OI(neuron),1))) = [];
%     end
% end

clear neuronIndex.OIindexTEMP
clearvars -except neuronIndex parameters_OI All
