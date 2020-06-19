function [psthBinsValue, PSTH, TrialDistr] = PSTH_event_SD_Cor_OmOnly(All, neuronIndex, parameters_OI, neuron, binSize, event, InfoDist, varargin)
%[psthBinsValue, PSTH, ID] = PSTH_event_SD(All, neuronIndex, parameters_OI, neuron, binSize, event,varargin)
% Function calculates the PSTHs, information distance measures, PSTH peaks
% and does bootstrap randomization
%
% V1.0
% h.terra@vu.nl
% May 16th 2019 
%
%

% Define time before and after events for PSTH calculation
TrialStart_tBef = 2;
Treshold_tBef = 2;
Response_tBef = 2;
TrialStart_tAft = 7;
Treshold_tAft = 3;
Response_tAft = 4;

  
    if strcmp(event,'All')
        eventCor = 'corCue';
        eventInc = 'incCue';
        eventOm = 'omCue';
        eventCorResp = 'corResp';
        eventIncResp = 'incResp';
        eventOmResp = 'omResp';
        eventPremResp = 'premResp';
        eventPremLatency = 'preLatency';
    elseif strcmp(event,'SD2')
        eventCor = 'corSD2Cue';
        eventInc = 'incSD2Cue';
        eventOm = 'omSD2Cue';
        eventCorResp = 'corSD2Resp';
        eventIncResp = 'incSD2Resp';
        eventOmResp = 'omSD2Resp';
        eventPremResp = 'premSD2Resp';
        eventPremLatency = 'preLatencySD2';
    elseif strcmp(event,'SD5')
        eventCor = 'corSD5Cue';
        eventInc = 'incSD5Cue';
        eventOm = 'omSD5Cue';
        eventCorResp = 'corSD5Resp';
        eventIncResp = 'incSD5Resp';
        eventOmResp = 'omSD5Resp';
        eventPremResp = 'premSD5Resp';
        eventPremLatency = 'preLatencySD5';
    elseif strcmp(event,'SD10')
        eventCor = 'corSD10Cue';
        eventInc = 'incSD10Cue';
        eventOm = 'omSD10Cue';
        eventCorResp = 'corSD10Resp';
        eventIncResp = 'incSD10Resp';
        eventOmResp = 'omSD10Resp';
        eventPremResp = 'premSD10Resp';
        eventPremLatency = 'preLatencySD10';
    end
    
    
    eventVideoBehavCor = eventVideo(All, eventCor, neuron, 0, 5,'PSTH_EVENT');
    eventVideoBehavInc = eventVideo(All, eventInc, neuron, 0, 5,'PSTH_EVENT');
    eventVideoBehavOm = eventVideo(All, eventOm, neuron, 0, 5,'PSTH_EVENT');
    eventVideoBehavPrem = eventVideo(All, eventPremResp, neuron, 0, 5,'PSTH_EVENT');
    eventTimes.Cor.trial_start = All{22, 1}{neuron, 1}.(eventCor)-5;
    eventTimes.Inc.trial_start = All{22, 1}{neuron, 1}.(eventInc)-5;
    eventTimes.Om.trial_start = All{22, 1}{neuron, 1}.(eventOm)-5;
    eventTimes.Prem.trial_start = All{22, 1}{neuron, 1}.(eventPremResp)-All{22, 1}{neuron, 1}.(eventPremLatency);
    eventTimes.Cor.wait_start = eventTimes.Cor.trial_start+eventVideoBehavCor.treshold';
    eventTimes.Inc.wait_start = eventTimes.Inc.trial_start+eventVideoBehavInc.treshold';
    eventTimes.Om.wait_start = eventTimes.Om.trial_start+eventVideoBehavOm.treshold';
    try
        eventTimes.Prem.wait_start = eventTimes.Prem.trial_start+eventVideoBehavPrem.treshold';
    catch
        eventTimes.Prem.wait_start = [];
    end
    eventTimes.Cor.resp = All{22, 1}{neuron, 1}.(eventCorResp);
    eventTimes.Inc.resp = All{22, 1}{neuron, 1}.(eventIncResp);
    eventTimes.Om.resp = All{22, 1}{neuron, 1}.(eventOmResp);
    eventTimes.Prem.resp = All{22, 1}{neuron, 1}.(eventPremResp);
    
    
    %[eventTimes_afterCorrect, eventTimes_afterError] = eventFilter_previousType(All, neuron);
    %eventTimes.Cor.cue = eventTimes_afterCorrect;
    
    
    % Get binned values per trial
    psthBinsValue.Cor.trial_start = [];
    psthBinsValue.Inc.trial_start = [];
    psthBinsValue.Om.trial_start = [];
    psthBinsValue.Prem.trial_start = [];
    psthBinsValue.Cor.wait_start = [];
    psthBinsValue.Inc.wait_start = [];
    psthBinsValue.Om.wait_start = [];
    psthBinsValue.Prem.wait_start = [];
    psthBinsValue.Cor.resp = [];
    psthBinsValue.Inc.resp = [];
    psthBinsValue.Om.resp = [];
    psthBinsValue.Prem.resp = [];
    
    % Find NaN values in eventTimes due to not crossing cue orientation
    % line, or video analysis not being trustworthy
    TrialDistr.CorTresholdExcluded = find(isnan(eventVideoBehavCor.treshold));
    TrialDistr.CorTresholdIncluded = find(~isnan(eventVideoBehavCor.treshold));
    TrialDistr.IncTresholdExcluded = find(isnan(eventVideoBehavInc.treshold));
    TrialDistr.IncTresholdIncluded = find(~isnan(eventVideoBehavInc.treshold));
    TrialDistr.OmTresholdExcluded = find(isnan(eventVideoBehavOm.treshold));
    TrialDistr.OmTresholdIncluded = find(~isnan(eventVideoBehavOm.treshold));
    TrialDistr.PremTresholdExcluded = find(isnan(eventVideoBehavPrem.treshold));
    TrialDistr.PremTresholdIncluded = find(~isnan(eventVideoBehavPrem.treshold));
    TrialDistr.CorTreshold = eventVideoBehavCor.treshold(TrialDistr.CorTresholdIncluded);
    TrialDistr.IncTreshold = eventVideoBehavInc.treshold(TrialDistr.IncTresholdIncluded);
    TrialDistr.OmTreshold = eventVideoBehavOm.treshold(TrialDistr.OmTresholdIncluded);
    TrialDistr.PremTreshold = eventVideoBehavPrem.treshold(TrialDistr.PremTresholdIncluded);
        
    psthBinsValue.Cor.trial_start = binData(eventTimes.Cor.trial_start(TrialDistr.CorTresholdIncluded), All{21, 1}{neuron, 1}, TrialStart_tBef, TrialStart_tAft, binSize);
    psthBinsValue.Inc.trial_start = binData(eventTimes.Inc.trial_start(TrialDistr.IncTresholdIncluded), All{21, 1}{neuron, 1}, TrialStart_tBef, TrialStart_tAft, binSize);
    psthBinsValue.Om.trial_start = binData(eventTimes.Om.trial_start(TrialDistr.OmTresholdIncluded), All{21, 1}{neuron, 1}, TrialStart_tBef, TrialStart_tAft, binSize);
    psthBinsValue.Prem.trial_start = binData(eventTimes.Prem.trial_start(TrialDistr.PremTresholdIncluded), All{21, 1}{neuron, 1}, TrialStart_tBef, TrialStart_tAft, binSize);
    psthBinsValue.Cor.wait_start = binData(eventTimes.Cor.wait_start(TrialDistr.CorTresholdIncluded), All{21, 1}{neuron, 1}, Treshold_tBef, Treshold_tAft, binSize);
    psthBinsValue.Inc.wait_start = binData(eventTimes.Inc.wait_start(TrialDistr.IncTresholdIncluded), All{21, 1}{neuron, 1}, Treshold_tBef, Treshold_tAft, binSize);
    psthBinsValue.Om.wait_start = binData(eventTimes.Om.wait_start(TrialDistr.OmTresholdIncluded), All{21, 1}{neuron, 1}, Treshold_tBef, Treshold_tAft, binSize);
    psthBinsValue.Prem.wait_start = binData(eventTimes.Prem.wait_start(TrialDistr.PremTresholdIncluded), All{21, 1}{neuron, 1}, Treshold_tBef, Treshold_tAft, binSize);
    psthBinsValue.Cor.resp = binData(eventTimes.Cor.resp(TrialDistr.CorTresholdIncluded), All{21, 1}{neuron, 1}, Response_tBef, Response_tAft, binSize);
    psthBinsValue.Inc.resp = binData(eventTimes.Inc.resp(TrialDistr.IncTresholdIncluded), All{21, 1}{neuron, 1}, Response_tBef, Response_tAft, binSize);
    psthBinsValue.Om.resp = binData(eventTimes.Om.resp(TrialDistr.OmTresholdIncluded), All{21, 1}{neuron, 1}, Response_tBef, Response_tAft, binSize);
    psthBinsValue.Prem.resp = binData(eventTimes.Prem.resp(TrialDistr.PremTresholdIncluded), All{21, 1}{neuron, 1}, Response_tBef, Response_tAft, binSize);
    
    %% Remove trials where animal went over treshold later than three seconds, passed the treshold again after first pass or had a lateral distance covered lower than 0.03 after treshold crossing
%     switch trialFilter
%         case 'yes'    
%             PSTH.Cor_ThreeSecsIx = find(eventVideoBehavCor.treshold <= 3 & eventVideoBehavCor.tresholdSecondPass == 0 & eventVideoBehavCor.tresholdLateralDist >= 0.03);
%             PSTH.Inc_ThreeSecsIx = find(eventVideoBehavInc.treshold <= 3 & eventVideoBehavInc.tresholdSecondPass == 0 & eventVideoBehavInc.tresholdLateralDist >= 0.03);
%             PSTH.Om_ThreeSecsIx = find(eventVideoBehavOm.treshold <= 3 & eventVideoBehavOm.tresholdSecondPass == 0 & eventVideoBehavOm.tresholdLateralDist >= 0.03);
%             PSTH.Prem_ThreeSecsIx = find(eventVideoBehavPrem.treshold <= 3 & eventVideoBehavPrem.tresholdSecondPass == 0 & eventVideoBehavPrem.tresholdLateralDist >= 0.03);
%             PSTH.Cor_ThreeSecsIxRemoved = find(eventVideoBehavCor.treshold > 3 | eventVideoBehavCor.tresholdSecondPass == 1 | eventVideoBehavCor.tresholdLateralDist < 0.03);
%             PSTH.Inc_ThreeSecsIxRemoved = find(eventVideoBehavInc.treshold > 3 | eventVideoBehavInc.tresholdSecondPass == 1 | eventVideoBehavInc.tresholdLateralDist < 0.03);
%             PSTH.Om_ThreeSecsIxRemoved = find(eventVideoBehavOm.treshold > 3 | eventVideoBehavOm.tresholdSecondPass == 1 | eventVideoBehavOm.tresholdLateralDist < 0.03);
%             PSTH.Prem_ThreeSecsIxRemoved = find(eventVideoBehavPrem.treshold > 3 | eventVideoBehavPrem.tresholdSecondPass == 1 | eventVideoBehavPrem.tresholdLateralDist < 0.03);
% 
%             %%
% 
%             psthBinsValue.Cor.trial_start = psthBinsValue.Cor.trial_start(PSTH.Cor_ThreeSecsIx,:);
%             psthBinsValue.Inc.trial_start = psthBinsValue.Inc.trial_start(PSTH.Inc_ThreeSecsIx,:);
%             psthBinsValue.Om.trial_start = psthBinsValue.Om.trial_start(PSTH.Om_ThreeSecsIx,:);
%             psthBinsValue.Prem.trial_start = psthBinsValue.Prem.trial_start(PSTH.Prem_ThreeSecsIx,:);
%             psthBinsValue.Cor.wait_start = psthBinsValue.Cor.wait_start(PSTH.Cor_ThreeSecsIx,:);
%             psthBinsValue.Inc.wait_start = psthBinsValue.Inc.wait_start(PSTH.Inc_ThreeSecsIx,:);
%             psthBinsValue.Om.wait_start = psthBinsValue.Om.wait_start(PSTH.Om_ThreeSecsIx,:);
%             psthBinsValue.Prem.wait_start = psthBinsValue.Prem.wait_start(PSTH.Prem_ThreeSecsIx,:);
%             psthBinsValue.Cor.resp = psthBinsValue.Cor.resp(PSTH.Cor_ThreeSecsIx,:);
%             psthBinsValue.Inc.resp = psthBinsValue.Inc.resp(PSTH.Inc_ThreeSecsIx,:);
%             psthBinsValue.Om.resp = psthBinsValue.Om.resp(PSTH.Om_ThreeSecsIx,:);
%             psthBinsValue.Prem.resp = psthBinsValue.Prem.resp(PSTH.Prem_ThreeSecsIx,:);
%         case 'no'
%     end


    
    nans = numel(find(isnan(psthBinsValue.Cor.trial_start(:,1))));
    PSTH.Cor.trial_start = (nansum(psthBinsValue.Cor.trial_start,1)/(size(psthBinsValue.Cor.trial_start,1)-nans))*(1/binSize);
    nans = numel(find(isnan(psthBinsValue.Inc.trial_start(:,1))));
    PSTH.Inc.trial_start = (nansum(psthBinsValue.Inc.trial_start,1)/(size(psthBinsValue.Inc.trial_start,1)-nans))*(1/binSize);
    nans = numel(find(isnan(psthBinsValue.Om.trial_start(:,1))));
    PSTH.Om.trial_start = (nansum(psthBinsValue.Om.trial_start,1)/(size(psthBinsValue.Om.trial_start,1)-nans))*(1/binSize);
    nans = numel(find(isnan(psthBinsValue.Prem.trial_start(:,1))));
    PSTH.Prem.trial_start = (nansum(psthBinsValue.Prem.trial_start,1)/(size(psthBinsValue.Prem.trial_start,1)-nans))*(1/binSize);
    
    nans = numel(find(isnan(psthBinsValue.Cor.wait_start(:,1))));
    PSTH.Cor.wait_start = (nansum(psthBinsValue.Cor.wait_start,1)/(size(psthBinsValue.Cor.wait_start,1)-nans))*(1/binSize);
    nans = numel(find(isnan(psthBinsValue.Inc.wait_start(:,1))));
    PSTH.Inc.wait_start = (nansum(psthBinsValue.Inc.wait_start,1)/(size(psthBinsValue.Inc.wait_start,1)-nans))*(1/binSize);
    nans = numel(find(isnan(psthBinsValue.Om.wait_start(:,1))));
    PSTH.Om.wait_start = (nansum(psthBinsValue.Om.wait_start,1)/(size(psthBinsValue.Om.wait_start,1)-nans))*(1/binSize);
    nans = numel(find(isnan(psthBinsValue.Prem.wait_start(:,1))));
    PSTH.Prem.wait_start = (nansum(psthBinsValue.Prem.wait_start,1)/(size(psthBinsValue.Prem.wait_start,1)-nans))*(1/binSize);
    
    nans = numel(find(isnan(psthBinsValue.Cor.resp(:,1))));
    PSTH.Cor.resp = (nansum(psthBinsValue.Cor.resp,1)/(size(psthBinsValue.Cor.resp,1)-nans))*(1/binSize);
    nans = numel(find(isnan(psthBinsValue.Inc.resp(:,1))));
    PSTH.Inc.resp = (nansum(psthBinsValue.Inc.resp,1)/(size(psthBinsValue.Inc.resp,1)-nans))*(1/binSize);
    nans = numel(find(isnan(psthBinsValue.Om.resp(:,1))));
    PSTH.Om.resp = (nansum(psthBinsValue.Om.resp,1)/(size(psthBinsValue.Om.resp,1)-nans))*(1/binSize);
    nans = numel(find(isnan(psthBinsValue.Prem.resp(:,1))));
    PSTH.Prem.resp = (nansum(psthBinsValue.Prem.resp,1)/(size(psthBinsValue.Prem.resp,1)-nans))*(1/binSize);
    
    %% Make PSTH from randomly shuffled event times and calculate peak and information distance
    switch InfoDist
        case 'yes'
            %Make randomized data
            outcome = {'Cor', 'Inc', 'Om', 'Prem'};
            for n = 1:numel(outcome)

                timePeriod = {'TrialStart', 'Treshold','Response'};
                for m = 1:numel(timePeriod)
                    if m == 1
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 0;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 5;
                        eventType = 'trial_start';
                        tBef = 2;
                        tAft = 7;
                    elseif m == 2
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 1;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 3;
                        eventType = 'wait_start';
                        tBef = 2;
                        tAft = 5;
                    elseif m == 3
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 0;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 3;
                        eventType = 'resp';
                        tBef = 2;
                        tAft = 4;
                    end


                BinsPerSec = 1/binSize;
                FRmean = mean(PSTH.Cor.trial_start((TrialStart_tBef*(1/binSize)-1*(1/binSize)):TrialStart_tBef*(1/binSize)));
                PSTH.ID.(outcome{n}).(timePeriod{m}).ID_real = InformDistance(FRmean, PSTH.(outcome{n}).(eventType), tBef, PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID, PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID, binSize);
                maxPeak = max(smoothdata(PSTH.(outcome{n}).(eventType)((tBef*BinsPerSec)-(PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID*BinsPerSec):(tBef*BinsPerSec)+(PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID*BinsPerSec))-All{15,1}(neuron,2),'gaussian',3));
                minPeak = min(smoothdata(PSTH.(outcome{n}).(eventType)((tBef*BinsPerSec)-(PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID*BinsPerSec):(tBef*BinsPerSec)+(PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID*BinsPerSec))-All{15,1}(neuron,2),'gaussian',3));
                PSTH.ID.(outcome{n}).(timePeriod{m}).Peak_real = max(abs([minPeak maxPeak]));
                if PSTH.ID.(outcome{n}).(timePeriod{m}).Peak_real == abs(minPeak)
                    PSTH.ID.(outcome{n}).(timePeriod{m}).Peak_PosOrNeg = -1;
                elseif PSTH.ID.(outcome{n}).(timePeriod{m}).Peak_real == abs(maxPeak)
                    PSTH.ID.(outcome{n}).(timePeriod{m}).Peak_PosOrNeg = 1;
                end

                % Do random analysis once and use that for multiple timePeriods.
                % Take over 1000 random timepoints ('events')
                bootstraps = 1000;
                PSTH.ID.(outcome{n}).(timePeriod{m}).Peak_rand = zeros(1,bootstraps);
                PSTH.ID.(outcome{n}).(timePeriod{m}).ID_rand = zeros(1,bootstraps);
                for bStrap = 1:bootstraps
                    eventNr = size(psthBinsValue.(outcome{n}).(eventType),1);
                    randEventTimes = rand(1,eventNr)*All{22, 1}{neuron, 1}.CH11(end-1);
                    %randEventTimes = rand(1,size(psthBinsValue.(outcome{n}).(eventType),1))*All{22, 1}{neuron, 1}.CH11(end-1);
                    BinsValueRand = binData(randEventTimes, All{21, 1}{neuron, 1}, tBef, tAft, binSize);
                    PSTHRand = (nansum(BinsValueRand,1)/(size(BinsValueRand,1)-nans))*(1/binSize);
                    maxPeak = max(smoothdata(PSTHRand((tBef*BinsPerSec)-(PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID*BinsPerSec):(tBef*BinsPerSec)+(PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID*BinsPerSec))-All{15,1}(neuron,2),'gaussian',3));
                    minPeak = min(smoothdata(PSTHRand((tBef*BinsPerSec)-(PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID*BinsPerSec):(tBef*BinsPerSec)+(PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID*BinsPerSec))-All{15,1}(neuron,2),'gaussian',3));
                    PSTH.ID.(outcome{n}).(timePeriod{m}).Peak_rand(bStrap) = max(abs([minPeak maxPeak]));
                    PSTH.ID.(outcome{n}).(timePeriod{m}).ID_rand(bStrap) = InformDistance(FRmean, PSTHRand, tBef, PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID, PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID, binSize);
                end

                PSTH.ID.(outcome{n}).(timePeriod{m}).ID_norm = PSTH.ID.(outcome{n}).(timePeriod{m}).ID_real/abs(mean(PSTH.ID.(outcome{n}).(timePeriod{m}).ID_rand));
                end
            end
        case 'no'
    end

    

    %% Plot OI raster, bar graph, waveform
    %% Waveform shape with neuron in black, position in mPFC, ISI
    %% PSTH TS, TH, Cue
    %% raster TS, TH, Cue
    %% 
    
    if nargin
    for iEvent=1:size(varargin,2)
        switch varargin{iEvent}
            case 'Plot'
                %% Get OI information
                % Get stimulus parameters
                Stim.Nr = [];
                Stim.Dur = All{1, 1}{neuron, 1}(:,1);
                Stim.Power = All{1, 1}{neuron, 1}(:,2);
                for i = 1:size(All{1, 1}{neuron, 1}(:,1),1)
                    Stim.Nr = [Stim.Nr 500];
                end
                lightIxs = find(All{1,1}{neuron,1}(:,3) < 0.01 & All{1,1}{neuron,1}(:,15)>= parameters_OI.PWCor & All{1,1}{neuron,1}(:,2)>= parameters_OI.minLightInt & All{1,1}{neuron,1}(:,6)>= parameters_OI.ReliabilityTreshold);
                [~, stimBlocksDur] = max(All{1,1}{neuron,1}(lightIxs,1));
                stimBlock = lightIxs(stimBlocksDur);
                %stimBlock = 1;
                if stimBlock == 1
                    range = stimBlock:Stim.Nr(stimBlock);
                else
                    range = 1+sum(Stim.Nr(1:stimBlock-1)):sum(Stim.Nr(1:stimBlock));
                end
                

                spikes = All{21, 1}{neuron,1};
                eventArrayOI = All{22, 1}{neuron,1}.CH16(range);
                if isempty(eventArrayOI)
                    eventArrayOI = 100;
                end

                %[spt_test, spt_baseline, FSLatency, jitter, reliability, spt_spikeIx, baselineTimeArray] = binMakerSALT2(All, eventArrayOI, spikes);

                
                figure
                
                
                    dim = [0.7 0.15 0.3 0.3];
                    content = {sprintf('Avg Spikes s^{-1} = %g',All{15,1}(neuron,2)), sprintf('ID %g', All{3, 1}(neuron,1)), sprintf('L-ratio %g', All{3, 1}(neuron,2)), sprintf('ISI1.5 %g',...
                    All{3, 1}(neuron,7)), sprintf('session dur (hour) %g', All{22, 1}{neuron, 1}.corResp(end)/3600), sprintf('shank %g', All{4, 1}(neuron,1)),...
                    sprintf('depth %g', All{11, 1}(neuron,1)), sprintf('animal %s', All{8, 1}{neuron,1}), sprintf('buddy %g', All{25, 1}(neuron,1)), sprintf('SALT_P %g', min(All{1, 1}{neuron,1}(:,3)))};
                    annotation('textbox',dim,'String',content,'FitBoxToText','on');
                    plotedit on
                if numel(find(ismember(neuronIndex.OIindex,neuron))) % Plot OI waveform if OI neuron is indexed
                    tBefOI = 0.04;
                    tAftOI = 0.04;
                    subplot(3,4,1);
                    for i=1:numel(eventArrayOI)
                        spnum=find(spikes>eventArrayOI(i)-tBefOI & spikes < eventArrayOI(i)+tAftOI);
                        scatter(spikes(spnum)-eventArrayOI(i),repmat(i,numel(spnum),1)',1,'k','filled');
                        hold on
                    end
                    axis([-tBefOI tAftOI 0 numel(eventArrayOI)])
                    title(sprintf('Light activation, %d ms, %d power', Stim.Dur(stimBlock), Stim.Power(stimBlock)));
                    ylabel('Trial')
                    xlabel('Time (ms)')

                    % get bins for binplot
                    binSizeOI = 0.001;
                    segmentSize = 80;

                    % Pre-allocate test and baseline matrix that are input for SALT test.
                    bins = zeros(length(eventArrayOI),segmentSize);

                    n=1;
                    for i = 1:length(eventArrayOI)
                        bins(i,:) = histcounts(spikes,eventArrayOI(i)-tBefOI:binSizeOI:eventArrayOI(i)+tAftOI);
                         %j = linspace(1,segmentSize,segmentSize);
%                         for m = 1:length(j)
%                             x = find(spikes>=((eventArrayOI(i)-tBefOI)+((j(m)-1)*binSizeOI)) & spikes<(eventArrayOI(i)-tBefOI)+((j(m)-1)*binSizeOI)+binSizeOI);
                            
%                             if ~isempty(x)
%                                 bins(i,j(m))=1;
%                                 n = n+1;
%                                 break
%                             end
%                         end
                    end
                    bins = sum(bins,1)./500;

                    subplot(3,4,2);
                    %histogram(bins,numel(-tBefOI:binSizeOI:tAftOI))
                    bar(linspace(-40,40,80),bins,'FaceColor','black');
                    xlim([-40 40]);
                    xlabel('Time from light onset (ms)')
                    ylabel('Avg. Spikes bin^{-1}')
                end
                
                subplot(3,4,3);
                hold on
                xt = linspace(0,(60)/(30003.0003/1000),61);
                shadedErrorBar(xt,All{13, 1}(neuron,45:105),All{14, 1}(neuron,45:105),{'k'});
                if numel(find(ismember(neuronIndex.OIindex,neuron))) % Plot OI waveform if OI neuron is indexed
                    try
                        avgWFOICor = All{16, 1}{neuron,stimBlock}(1,45:105);
                    catch
                        avgWFOICor = All{13, 1}(neuron,45:105);
                    end
                    plot(xt, avgWFOICor,'b','linewidth',2)
                end
                ylabel('Amplitude (uV)')
                xlabel('time (ms)')
                title(sprintf('neuron %d', neuron));
                
                subplot(3,4,4)
                spikeIntervals = diff(All{21,1}{neuron,1}(1:All{3, 1}(neuron,3)))*1000;
                binSize_hist = 0.5;                                            % 1 ms bins
                x = [0:binSize_hist:50];
                histogram(spikeIntervals,x,'FaceColor','black')
                ylabel('Count')
                xlabel('ISI (ms)')
                xlim([0 50])
                
                % Get timeaxis for PSTH
                time_axis_Treshold = round(linspace(-Treshold_tBef+(binSize),Treshold_tAft,((Treshold_tBef+Treshold_tAft)/binSize)),1);
                time_axis_TrialStart = round(linspace(-TrialStart_tBef+(binSize),TrialStart_tAft,((TrialStart_tBef+TrialStart_tAft)/binSize)),1);
                time_axis_Response = round(linspace(-Response_tBef+(binSize),Response_tAft,((Response_tBef+Response_tAft)/binSize)),1);
                
                subplot(3,4,5)
                hold on
                SEM = nanstd(psthBinsValue.Cor.trial_start,1)/sqrt(size(psthBinsValue.Cor.trial_start,1)); 
                shadedErrorBar(time_axis_TrialStart,PSTH.Cor.trial_start,SEM,{'k'});
%                 SEM = nanstd(psthBinsValue.Inc.trial_start,1)/sqrt(size(psthBinsValue.Inc.trial_start,1)); 
%                 shadedErrorBar(time_axis_TrialStart,PSTH.Inc.trial_start,SEM,{'r'});
                SEM = nanstd(psthBinsValue.Om.trial_start,1)/sqrt(size(psthBinsValue.Om.trial_start,1)); 
                shadedErrorBar(time_axis_TrialStart,PSTH.Om.trial_start,SEM,{'b'});
%                 SEM = nanstd(psthBinsValue.Prem.trial_start,1)/sqrt(size(psthBinsValue.Prem.trial_start,1)); 
%                 shadedErrorBar(time_axis_TrialStart,PSTH.Prem.trial_start,SEM,{'c'});
%                 axis([-TrialStart_tBef TrialStart_tAft 0 inf])
                
                subplot(3,4,6)
                hold on
                SEM = nanstd(psthBinsValue.Cor.wait_start,1)/sqrt(size(psthBinsValue.Cor.wait_start,1)); 
                shadedErrorBar(time_axis_Treshold,PSTH.Cor.wait_start,SEM,{'k'})
%                 SEM = nanstd(psthBinsValue.Inc.wait_start,1)/sqrt(size(psthBinsValue.Inc.wait_start,1)); 
%                 shadedErrorBar(time_axis_Treshold,PSTH.Inc.wait_start,SEM,{'r'});
                SEM = nanstd(psthBinsValue.Om.wait_start,1)/sqrt(size(psthBinsValue.Om.wait_start,1)); 
                shadedErrorBar(time_axis_Treshold,PSTH.Om.wait_start,SEM,{'b'});
%                 SEM = nanstd(psthBinsValue.Prem.wait_start,1)/sqrt(size(psthBinsValue.Prem.wait_start,1)); 
%                 shadedErrorBar(time_axis_Treshold,PSTH.Prem.wait_start,SEM,{'c'});
%                 axis([-Treshold_tBef Treshold_tAft 0 inf])

                subplot(3,4,7)
                SEM = nanstd(psthBinsValue.Cor.resp,1)/sqrt(size(psthBinsValue.Cor.resp,1)); 
                shadedErrorBar(time_axis_Response,PSTH.Cor.resp,SEM,{'k'});
%                 hold on
%                 SEM = nanstd(psthBinsValue.Inc.resp,1)/sqrt(size(psthBinsValue.Inc.resp,1)); 
%                 shadedErrorBar(time_axis_Response,PSTH.Inc.resp,SEM,{'r'});
%                 hold on
                SEM = nanstd(psthBinsValue.Om.resp,1)/sqrt(size(psthBinsValue.Om.resp,1)); 
                shadedErrorBar(time_axis_Response,PSTH.Om.resp,SEM,{'b'});
%                 SEM = nanstd(psthBinsValue.Prem.resp,1)/sqrt(size(psthBinsValue.Prem.resp,1)); 
%                 shadedErrorBar(time_axis_Response,PSTH.Prem.resp,SEM,{'c'});
%                 axis([-Response_tBef Response_tAft 0 inf])
                
                dotSize = 0.3;
                subplot(3,4,9);
                hold on
                n = 1;
                for i=1:numel(TrialDistr.CorTresholdIncluded)
                    spnum=find(All{21, 1}{neuron, 1}>eventTimes.Cor.trial_start(TrialDistr.CorTresholdIncluded(i))-TrialStart_tBef & All{21, 1}{neuron, 1} < eventTimes.Cor.trial_start(TrialDistr.CorTresholdIncluded(i))+TrialStart_tAft);
                    scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Cor.trial_start(TrialDistr.CorTresholdIncluded(i)),repmat(n,numel(spnum),1)',dotSize,'k','filled');
                    scatter(eventVideoBehavCor.treshold(i),repmat(n,numel(1),1)',dotSize*2,'m','filled');
                    n = n+1;
                end
% 
%                 for i=1:numel(TrialDistr.IncTresholdIncluded)
%                     spnum=find(All{21, 1}{neuron, 1}>eventTimes.Inc.trial_start(TrialDistr.IncTresholdIncluded(i))-TrialStart_tBef & All{21, 1}{neuron, 1} < eventTimes.Inc.trial_start(TrialDistr.IncTresholdIncluded(i))+TrialStart_tAft);
%                     scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Inc.trial_start(TrialDistr.IncTresholdIncluded(i)),repmat(n,numel(spnum),1)',0.5,'r','filled');
%                     n = n+1;
%                 end
% 
                for i=1:numel(TrialDistr.OmTresholdIncluded)
                    spnum=find(All{21, 1}{neuron, 1}>eventTimes.Om.trial_start(TrialDistr.OmTresholdIncluded(i))-TrialStart_tBef & All{21, 1}{neuron, 1} < eventTimes.Om.trial_start(TrialDistr.OmTresholdIncluded(i))+TrialStart_tAft);
                    scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Om.trial_start(TrialDistr.OmTresholdIncluded(i)),repmat(n,numel(spnum),1)',dotSize,'b','filled');
                    scatter(eventVideoBehavOm.treshold(i),repmat(n,numel(1),1)',dotSize*2,'m','filled');
                    n = n+1;
                end
%                 for i=1:numel(TrialDistr.PremTresholdIncluded)
%                     spnum=find(All{21, 1}{neuron, 1}>eventTimes.Prem.trial_start(TrialDistr.PremTresholdIncluded(i))-TrialStart_tBef & All{21, 1}{neuron, 1} < eventTimes.Prem.trial_start(TrialDistr.PremTresholdIncluded(i))+TrialStart_tAft);
%                     scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Prem.trial_start(TrialDistr.PremTresholdIncluded(i)),repmat(n,numel(spnum),1)',0.5,'c','filled');
%                     n = n+1;
%                 end
                axis([-TrialStart_tBef TrialStart_tAft 1 numel(TrialDistr.CorTresholdIncluded)+numel(TrialDistr.OmTresholdIncluded)])
                ylabel('Trial')
                xlabel('Time start trial (sec)')
                subplot(3,4,10);
                hold on
                n = 1;
                for i=1:numel(TrialDistr.CorTresholdIncluded)
                    spnum=find(All{21, 1}{neuron, 1}>eventTimes.Cor.wait_start(TrialDistr.CorTresholdIncluded(i))-Treshold_tBef & All{21, 1}{neuron, 1} < eventTimes.Cor.wait_start(TrialDistr.CorTresholdIncluded(i))+Treshold_tAft);
                    scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Cor.wait_start(TrialDistr.CorTresholdIncluded(i)),repmat(n,numel(spnum),1)',dotSize,'k','filled');
                    scatter(eventVideoBehavCor.treshold(i),repmat(n,numel(1),1)',dotSize*2,'m','filled');
                    n = n+1;
                end
% 
%                 for i=1:numel(TrialDistr.IncTresholdIncluded)
%                     spnum=find(All{21, 1}{neuron, 1}>eventTimes.Inc.wait_start(TrialDistr.IncTresholdIncluded(i))-Treshold_tBef & All{21, 1}{neuron, 1} < eventTimes.Inc.wait_start(TrialDistr.IncTresholdIncluded(i))+Treshold_tAft);
%                     scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Inc.wait_start(TrialDistr.IncTresholdIncluded(i)),repmat(n,numel(spnum),1)',0.5,'r','filled');
%                     n = n+1;
%                 end
                for i=1:numel(TrialDistr.OmTresholdIncluded)
                    spnum=find(All{21, 1}{neuron, 1}>eventTimes.Om.wait_start(TrialDistr.OmTresholdIncluded(i))-Treshold_tBef & All{21, 1}{neuron, 1} < eventTimes.Om.wait_start(TrialDistr.OmTresholdIncluded(i))+Treshold_tAft);
                    scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Om.wait_start(TrialDistr.OmTresholdIncluded(i)),repmat(n,numel(spnum),1)',dotSize,'b','filled');
                    n = n+1;
                end
%                 for i=1:numel(TrialDistr.PremTresholdIncluded)
%                     spnum=find(All{21, 1}{neuron, 1}>eventTimes.Prem.wait_start(TrialDistr.PremTresholdIncluded(i))-Treshold_tBef & All{21, 1}{neuron, 1} < eventTimes.Prem.wait_start(TrialDistr.PremTresholdIncluded(i))+Treshold_tAft);
%                     scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Prem.wait_start(TrialDistr.PremTresholdIncluded(i)),repmat(n,numel(spnum),1)',0.5,'c','filled');
%                     n = n+1;
%                 end
                axis([-Treshold_tBef Treshold_tAft 1 numel(TrialDistr.CorTresholdIncluded)+numel(TrialDistr.OmTresholdIncluded)])
                ylabel('Trial')
                xlabel('Time cue orientation (sec)')

                subplot(3,4,11);
                hold on
                n = 1;
                for i=1:numel(TrialDistr.CorTresholdIncluded)
                    spnum=find(All{21, 1}{neuron, 1}>eventTimes.Cor.resp(TrialDistr.CorTresholdIncluded(i))-Response_tBef & All{21, 1}{neuron, 1} < eventTimes.Cor.resp(TrialDistr.CorTresholdIncluded(i))+Response_tAft);
                    scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Cor.resp(TrialDistr.CorTresholdIncluded(i)),repmat(n,numel(spnum),1)',dotSize,'k','filled');
                    scatter(eventVideoBehavCor.treshold(i),repmat(n,numel(1),1)',dotSize*2,'m','filled');
                    n = n+1;
                end
% 
%                 for i=1:numel(TrialDistr.IncTresholdIncluded)
%                     spnum=find(All{21, 1}{neuron, 1}>eventTimes.Inc.resp(TrialDistr.IncTresholdIncluded(i))-Response_tBef & All{21, 1}{neuron, 1} < eventTimes.Inc.resp(TrialDistr.IncTresholdIncluded(i))+Response_tAft);
%                     scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Inc.resp(TrialDistr.IncTresholdIncluded(i)),repmat(n,numel(spnum),1)',0.5,'r','filled');
%                     n = n+1;
%                 end
% 
                for i=1:numel(TrialDistr.OmTresholdIncluded)
                    spnum=find(All{21, 1}{neuron, 1}>eventTimes.Om.resp(TrialDistr.OmTresholdIncluded(i))-Response_tBef & All{21, 1}{neuron, 1} < eventTimes.Om.resp(TrialDistr.OmTresholdIncluded(i))+Response_tAft);
                    scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Om.resp(TrialDistr.OmTresholdIncluded(i)),repmat(n,numel(spnum),1)',dotSize,'b','filled');
                    n = n+1;
                end
%                 for i=1:numel(TrialDistr.PremTresholdIncluded)
%                     spnum=find(All{21, 1}{neuron, 1}>eventTimes.Prem.resp(TrialDistr.PremTresholdIncluded(i))-Response_tBef & All{21, 1}{neuron, 1} < eventTimes.Prem.resp(TrialDistr.PremTresholdIncluded(i))+Response_tAft);
%                     scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Prem.resp(TrialDistr.PremTresholdIncluded(i)),repmat(n,numel(spnum),1)',0.5,'c','filled');
%                     n = n+1;
%                 end
                axis([-Response_tBef Response_tAft 1 numel(TrialDistr.CorTresholdIncluded)+numel(TrialDistr.OmTresholdIncluded)])
                ylabel('Trial')
                xlabel('Time response (sec)')
                
                %formatSpec = 'VARSD_all_neuron%d.fig';
                %str = sprintf(formatSpec, neuron);
                %savefig(str)
                %print(str,'-deps','-painters')
                %movefile (str, 'var_SD_all');
                %close
        end
    end
    end