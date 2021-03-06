function [psthBinsValue, PSTH, TrialDistr] = PSTH_event_ITI125(All, neuronIndex, parameters_OI, neuron, binSize, event, InfoDist,varargin)

% Figure 6A and 6B


[Correct, Incorrect, Omission, TrialDistr.Premature, Perserverative, nrTrialsStage, Perf]  = CC_ephys_var_ITI_SD_TTL_Master(10,1,All{5, 1}{neuron, 2});   
if strcmp(event,'All')
    [Correct, Incorrect, Omission, Premature, Perserverative, nrTrialsStage, Perf]  = CC_ephys_var_ITI_SD_TTL_Master(10,1,All{5, 1}{neuron, 2});
    eventCor = 'corCue';
    eventInc = 'incCue';
    eventOm = 'omCue';
    eventCorResp = 'corResp';
    eventIncResp = 'incResp';
    eventOmResp = 'omResp';
    eventPremResp = 'premResp';
    eventTimes.Cor.cue = All{22, 1}{neuron, 1}.(eventCor);
    eventTimes.Inc.cue = All{22, 1}{neuron, 1}.(eventInc);
    eventTimes.Om.cue = All{22, 1}{neuron, 1}.(eventOm);

    eventTimes.Cor.trial_start = zeros(1,Correct.TotalNr);
    eventTimes.Inc.trial_start = zeros(1,Incorrect.TotalNr);
    eventTimes.Om.trial_start = zeros(1,Omission.TotalNr);
    eventTimes.Cor.trial_start(Correct.TrialPosITI50) = All{22, 1}{neuron, 1}.(eventCor)(Correct.TrialPosITI50)-5.0;
    eventTimes.Inc.trial_start(Incorrect.TrialPosITI50) = All{22, 1}{neuron, 1}.(eventInc)(Incorrect.TrialPosITI50)-5.0;
    eventTimes.Om.trial_start(Omission.TrialPosITI50) = All{22, 1}{neuron, 1}.(eventOm)(Omission.TrialPosITI50)-5.0;
    eventTimes.Cor.trial_start(Correct.TrialPosITI75) = All{22, 1}{neuron, 1}.(eventCor)(Correct.TrialPosITI75)-7.5;
    eventTimes.Inc.trial_start(Incorrect.TrialPosITI75) = All{22, 1}{neuron, 1}.(eventInc)(Incorrect.TrialPosITI75)-7.5;
    eventTimes.Om.trial_start(Omission.TrialPosITI75) = All{22, 1}{neuron, 1}.(eventOm)(Omission.TrialPosITI75)-7.5;
    eventTimes.Cor.trial_start(Correct.TrialPosITI125) = All{22, 1}{neuron, 1}.(eventCor)(Correct.TrialPosITI125)-12.5;
    eventTimes.Inc.trial_start(Incorrect.TrialPosITI125) = All{22, 1}{neuron, 1}.(eventInc)(Incorrect.TrialPosITI125)-12.5;
    eventTimes.Om.trial_start(Omission.TrialPosITI125) = All{22, 1}{neuron, 1}.(eventOm)(Omission.TrialPosITI125)-12.5;
    eventTimes.Prem.trial_start = All{22, 1}{neuron, 1}.(eventPremResp)-All{22, 1}{neuron, 1}.preLatency;

    eventVideoBehavCor.treshold = zeros(1,Correct.TotalNr);
    eventVideoBehavCor.tresholdSecondPass = zeros(1,Correct.TotalNr);
    eventVideoBehavCor.tresholdLateralDist = zeros(1,Correct.TotalNr);
    eventVideoBehavInc.treshold = zeros(1,Incorrect.TotalNr);
    eventVideoBehavInc.tresholdSecondPass = zeros(1,Incorrect.TotalNr);
    eventVideoBehavInc.tresholdLateralDist = zeros(1,Incorrect.TotalNr);
    eventVideoBehavOm.treshold = zeros(1,Omission.TotalNr);
    eventVideoBehavOm.tresholdSecondPass = zeros(1,Omission.TotalNr);
    eventVideoBehavOm.tresholdLateralDist = zeros(1,Omission.TotalNr);
    eventVideoBehavPrem.treshold = zeros(1,Premature.TotalNr);
    eventVideoBehavPrem.tresholdSecondPass = zeros(1,Premature.TotalNr);
    eventVideoBehavPrem.tresholdLateralDist = zeros(1,Premature.TotalNr);

    eventVideoBehavCorITI50 = eventVideo(All, 'corITI50Cue', neuron, 0, 5.0,'PSTH_EVENT');
    eventVideoBehavCor.treshold(Correct.TrialPosITI50) = eventVideoBehavCorITI50.treshold;
    eventVideoBehavCor.tresholdSecondPass(Correct.TrialPosITI50) = eventVideoBehavCorITI50.tresholdSecondPass;
    eventVideoBehavCor.tresholdLateralDist(Correct.TrialPosITI50) = eventVideoBehavCorITI50.tresholdLateralDist;
    eventVideoBehavCorITI75 = eventVideo(All, 'corITI75Cue', neuron, 0, 7.5,'PSTH_EVENT');
    eventVideoBehavCor.treshold(Correct.TrialPosITI75) = eventVideoBehavCorITI75.treshold;
    eventVideoBehavCor.tresholdSecondPass(Correct.TrialPosITI75) = eventVideoBehavCorITI75.tresholdSecondPass;
    eventVideoBehavCor.tresholdLateralDist(Correct.TrialPosITI75) = eventVideoBehavCorITI75.tresholdLateralDist;
    eventVideoBehavCorITI125 = eventVideo(All, 'corITI125Cue', neuron, 0, 12.5,'PSTH_EVENT');
    eventVideoBehavCor.treshold(Correct.TrialPosITI125) = eventVideoBehavCorITI125.treshold;
    eventVideoBehavCor.tresholdSecondPass(Correct.TrialPosITI125) = eventVideoBehavCorITI125.tresholdSecondPass;
    eventVideoBehavCor.tresholdLateralDist(Correct.TrialPosITI125) = eventVideoBehavCorITI125.tresholdLateralDist;

    eventVideoBehavIncITI50 = eventVideo(All, 'incITI50Cue', neuron, 0, 5.0,'PSTH_EVENT');
    eventVideoBehavInc.treshold(Incorrect.TrialPosITI50) = eventVideoBehavIncITI50.treshold;
    eventVideoBehavInc.tresholdSecondPass(Incorrect.TrialPosITI50) = eventVideoBehavIncITI50.tresholdSecondPass;
    eventVideoBehavInc.tresholdLateralDist(Incorrect.TrialPosITI50) = eventVideoBehavIncITI50.tresholdLateralDist;
    eventVideoBehavIncITI75 = eventVideo(All, 'incITI75Cue', neuron, 0, 7.5,'PSTH_EVENT');
    eventVideoBehavInc.treshold(Incorrect.TrialPosITI75) = eventVideoBehavIncITI75.treshold;
    eventVideoBehavInc.tresholdSecondPass(Incorrect.TrialPosITI75) = eventVideoBehavIncITI75.tresholdSecondPass;
    eventVideoBehavInc.tresholdLateralDist(Incorrect.TrialPosITI75) = eventVideoBehavIncITI75.tresholdLateralDist;
    eventVideoBehavIncITI125 = eventVideo(All, 'incITI125Cue', neuron, 0, 12.5,'PSTH_EVENT');
    eventVideoBehavInc.treshold(Incorrect.TrialPosITI125) = eventVideoBehavIncITI125.treshold;
    eventVideoBehavInc.tresholdSecondPass(Incorrect.TrialPosITI125) = eventVideoBehavIncITI125.tresholdSecondPass;
    eventVideoBehavInc.tresholdLateralDist(Incorrect.TrialPosITI125) = eventVideoBehavIncITI125.tresholdLateralDist;

    eventVideoBehavOmITI50 = eventVideo(All, 'omITI50Cue', neuron, 0, 5.0,'PSTH_EVENT');
    eventVideoBehavOm.treshold(Omission.TrialPosITI50) = eventVideoBehavOmITI50.treshold;
    eventVideoBehavOm.tresholdSecondPass(Omission.TrialPosITI50) = eventVideoBehavOmITI50.tresholdSecondPass;
    eventVideoBehavOm.tresholdLateralDist(Omission.TrialPosITI50) = eventVideoBehavOmITI50.tresholdLateralDist;
    eventVideoBehavOmITI75 = eventVideo(All, 'omITI75Cue', neuron, 0, 7.5,'PSTH_EVENT');
    eventVideoBehavOm.treshold(Omission.TrialPosITI75) = eventVideoBehavOmITI75.treshold;
    eventVideoBehavOm.tresholdSecondPass(Omission.TrialPosITI75) = eventVideoBehavOmITI75.tresholdSecondPass;
    eventVideoBehavOm.tresholdLateralDist(Omission.TrialPosITI75) = eventVideoBehavOmITI75.tresholdLateralDist;
    eventVideoBehavOmITI125 = eventVideo(All, 'omITI125Cue', neuron, 0, 12.5,'PSTH_EVENT');
    eventVideoBehavOm.treshold(Omission.TrialPosITI125) = eventVideoBehavOmITI125.treshold;
    eventVideoBehavOm.tresholdSecondPass(Omission.TrialPosITI125) = eventVideoBehavOmITI125.tresholdSecondPass;
    eventVideoBehavOm.tresholdLateralDist(Omission.TrialPosITI125) = eventVideoBehavOmITI125.tresholdLateralDist;

    try
        eventVideoBehavPremITI50 = eventVideo(All, 'premITI50Resp', neuron, 0, 5.0,'PSTH_EVENT');
        eventVideoBehavPrem.treshold(Premature.TrialPosITI50) = eventVideoBehavPremITI50.treshold;
        eventVideoBehavPrem.tresholdSecondPass(Premature.TrialPosITI50) = eventVideoBehavPremITI50.tresholdSecondPass;
        eventVideoBehavPrem.tresholdLateralDist(Premature.TrialPosITI50) = eventVideoBehavPremITI50.tresholdLateralDist;
    catch
        eventVideoBehavPrem.treshold = [];
        eventVideoBehavPrem.tresholdSecondPass = [];
        eventVideoBehavPrem.tresholdLateralDist = [];
    end
    eventVideoBehavPremITI75 = eventVideo(All, 'premITI75Resp', neuron, 0, 7.5,'PSTH_EVENT');
    eventVideoBehavPrem.treshold(Premature.TrialPosITI75) = eventVideoBehavPremITI75.treshold;
    eventVideoBehavPrem.tresholdSecondPass(Premature.TrialPosITI75) = eventVideoBehavPremITI75.tresholdSecondPass;
    eventVideoBehavPrem.tresholdLateralDist(Premature.TrialPosITI75) = eventVideoBehavPremITI75.tresholdLateralDist;
    eventVideoBehavPremITI125 = eventVideo(All, 'premITI125Resp', neuron, 0, 12.5,'PSTH_EVENT');
    eventVideoBehavPrem.treshold(Premature.TrialPosITI125) = eventVideoBehavPremITI125.treshold;
    eventVideoBehavPrem.tresholdSecondPass(Premature.TrialPosITI125) = eventVideoBehavPremITI125.tresholdSecondPass;
    eventVideoBehavPrem.tresholdLateralDist(Premature.TrialPosITI125) = eventVideoBehavPremITI125.tresholdLateralDist;

    Cue_tBef = 2;
    TrialStart_tBef = 2;
    Treshold_tBef = 3;
    Response_tBef = 3;
    Cue_tAft = 2;
    TrialStart_tAft = 7;
    Treshold_tAft = 5;
    Response_tAft = 4;
    
    eventTimes.Cor.wait_start = eventTimes.Cor.trial_start+eventVideoBehavCor.treshold;
    eventTimes.Inc.wait_start = eventTimes.Inc.trial_start+eventVideoBehavInc.treshold;
    eventTimes.Om.wait_start = eventTimes.Om.trial_start+eventVideoBehavOm.treshold;
    try
        eventTimes.Prem.wait_start = eventTimes.Prem.trial_start+eventVideoBehavPrem.treshold;

    catch
        eventTimes.Prem.wait_start = [];
    end
    eventTimes.Cor.cue = All{22, 1}{neuron, 1}.(eventCor);
    eventTimes.Cor.resp = All{22, 1}{neuron, 1}.(eventCorResp);
    eventTimes.Inc.cue = All{22, 1}{neuron, 1}.(eventInc);
    eventTimes.Inc.resp = All{22, 1}{neuron, 1}.(eventIncResp);
    eventTimes.Om.cue = All{22, 1}{neuron, 1}.(eventOm);
    eventTimes.Om.resp = All{22, 1}{neuron, 1}.(eventOm);
    try
        eventTimes.Prem.resp = All{22, 1}{neuron, 1}.(eventPremResp);
    catch
        eventTimes.Prem.resp = [];
    end
elseif strcmp(event,'ITI50')
    eventCor = 'corITI50Cue';
    eventInc = 'incITI50Cue';
    eventOm = 'omITI50Cue';
    eventCorResp = 'corITI50Resp';
    eventIncResp = 'incITI50Resp';
    eventOmResp = 'omITI50Resp';
    eventPremResp = 'premITI50Resp';
    eventTimes.Cor.trial_start = All{22, 1}{neuron, 1}.(eventCor)-5;
    eventTimes.Inc.trial_start = All{22, 1}{neuron, 1}.(eventInc)-5;
    eventTimes.Om.trial_start = All{22, 1}{neuron, 1}.(eventOm)-5;
    eventTimes.Prem.trial_start = All{22, 1}{neuron, 1}.(eventPremResp)-All{22, 1}{neuron, 1}.preLatencyITI5;
    eventVideoBehavCor = eventVideo(All, eventCor, neuron, 0, 5,'PSTH_EVENT');
    eventVideoBehavInc = eventVideo(All, eventInc, neuron, 0, 5,'PSTH_EVENT');
    eventVideoBehavOm = eventVideo(All, eventOm, neuron, 0, 5,'PSTH_EVENT');
    eventVideoBehavPrem = eventVideo(All, eventPremResp, neuron, 0, 5,'PSTH_EVENT');
    Cue_tBef = 2;
    TrialStart_tBef = 2;
    Treshold_tBef = 2;
    Response_tBef = 2;
    Cue_tAft = 2;
    TrialStart_tAft = 7;
    Treshold_tAft = 5;
    Response_tAft = 4;
    
        eventTimes.Cor.wait_start = eventTimes.Cor.trial_start+eventVideoBehavCor.treshold';
    eventTimes.Inc.wait_start = eventTimes.Inc.trial_start+eventVideoBehavInc.treshold';
    eventTimes.Om.wait_start = eventTimes.Om.trial_start+eventVideoBehavOm.treshold';
    try
        eventTimes.Prem.wait_start = eventTimes.Prem.trial_start+eventVideoBehavPrem.treshold';

    catch
        eventTimes.Prem.wait_start = [];
    end
    eventTimes.Cor.cue = All{22, 1}{neuron, 1}.(eventCor);
    eventTimes.Cor.resp = All{22, 1}{neuron, 1}.(eventCorResp);
    eventTimes.Inc.cue = All{22, 1}{neuron, 1}.(eventInc);
    eventTimes.Inc.resp = All{22, 1}{neuron, 1}.(eventIncResp);
    eventTimes.Om.cue = All{22, 1}{neuron, 1}.(eventOm);
    eventTimes.Om.resp = All{22, 1}{neuron, 1}.(eventOm);
    try
        eventTimes.Prem.resp = All{22, 1}{neuron, 1}.(eventPremResp);
    catch
        eventTimes.Prem.resp = [];
    end
    
elseif strcmp(event,'ITI75')
    eventCor = 'corITI75Cue';
    eventInc = 'incITI75Cue';
    eventOm = 'omITI75Cue';
    eventCorResp = 'corITI75Resp';
    eventIncResp = 'incITI75Resp';
    eventOmResp = 'omITI75Resp';
    eventPremResp = 'premITI75Resp';
    eventTimes.Cor.trial_start = All{22, 1}{neuron, 1}.(eventCor)-7.5;
    eventTimes.Inc.trial_start = All{22, 1}{neuron, 1}.(eventInc)-7.5;
    eventTimes.Om.trial_start = All{22, 1}{neuron, 1}.(eventOm)-7.5;
    eventTimes.Prem.trial_start = All{22, 1}{neuron, 1}.(eventPremResp)-All{22, 1}{neuron, 1}.preLatencyITI75;
    eventVideoBehavCor = eventVideo(All, eventCor, neuron, 0, 7.5,'PSTH_EVENT');
    eventVideoBehavInc = eventVideo(All, eventInc, neuron, 0, 7.5,'PSTH_EVENT');
    eventVideoBehavOm = eventVideo(All, eventOm, neuron, 0, 7.5,'PSTH_EVENT');
    eventVideoBehavPrem = eventVideo(All, eventPremResp, neuron, 0, 7.5,'PSTH_EVENT');
    Cue_tBef = 2;
    TrialStart_tBef = 2;
    Treshold_tBef = 2;
    Response_tBef = 2;
    Cue_tAft = 2;
    TrialStart_tAft = 10;
    Treshold_tAft = 8;
    Response_tAft = 4;
    
        eventTimes.Cor.wait_start = eventTimes.Cor.trial_start+eventVideoBehavCor.treshold';
    eventTimes.Inc.wait_start = eventTimes.Inc.trial_start+eventVideoBehavInc.treshold';
    eventTimes.Om.wait_start = eventTimes.Om.trial_start+eventVideoBehavOm.treshold';
    try
        eventTimes.Prem.wait_start = eventTimes.Prem.trial_start+eventVideoBehavPrem.treshold';

    catch
        eventTimes.Prem.wait_start = [];
    end
    eventTimes.Cor.cue = All{22, 1}{neuron, 1}.(eventCor);
    eventTimes.Cor.resp = All{22, 1}{neuron, 1}.(eventCorResp);
    eventTimes.Inc.cue = All{22, 1}{neuron, 1}.(eventInc);
    eventTimes.Inc.resp = All{22, 1}{neuron, 1}.(eventIncResp);
    eventTimes.Om.cue = All{22, 1}{neuron, 1}.(eventOm);
    eventTimes.Om.resp = All{22, 1}{neuron, 1}.(eventOm);
    try
        eventTimes.Prem.resp = All{22, 1}{neuron, 1}.(eventPremResp);
    catch
        eventTimes.Prem.resp = [];
    end
elseif strcmp(event,'ITI125')
    eventCor = 'corITI125Cue';
    eventInc = 'incITI125Cue';
    eventOm = 'omITI125Cue';
    eventCorResp = 'corITI125Resp';
    eventIncResp = 'incITI125Resp';
    eventOmResp = 'omITI125Resp';
    eventPremResp = 'premITI125Resp';
    eventTimes.Cor.trial_start = All{22, 1}{neuron, 1}.(eventCor)-12.5;
    eventTimes.Inc.trial_start = All{22, 1}{neuron, 1}.(eventInc)-12.5;
    eventTimes.Om.trial_start = All{22, 1}{neuron, 1}.(eventOm)-12.5;
    eventTimes.Prem.trial_start = All{22, 1}{neuron, 1}.(eventPremResp)-All{22, 1}{neuron, 1}.preLatencyITI125;
    eventVideoBehavCor = eventVideo(All, eventCor, neuron, 0, 12.5,'PSTH_EVENT');
    eventVideoBehavInc = eventVideo(All, eventInc, neuron, 0, 12.5,'PSTH_EVENT');
    eventVideoBehavOm = eventVideo(All, eventOm, neuron, 0, 12.5,'PSTH_EVENT');
    eventVideoBehavPrem = eventVideo(All, eventPremResp, neuron, 0, 12.5,'PSTH_EVENT');
    Cue_tBef = 2;
    TrialStart_tBef = 2;
    Treshold_tBef = 2;
    Response_tBef = 2;
    Cue_tAft = 2;
    TrialStart_tAft = 15;
    Treshold_tAft = 13;
    Response_tAft = 4;
    
        eventTimes.Cor.wait_start = eventTimes.Cor.trial_start+eventVideoBehavCor.treshold';
    eventTimes.Inc.wait_start = eventTimes.Inc.trial_start+eventVideoBehavInc.treshold';
    eventTimes.Om.wait_start = eventTimes.Om.trial_start+eventVideoBehavOm.treshold';
    try
        eventTimes.Prem.wait_start = eventTimes.Prem.trial_start+eventVideoBehavPrem.treshold';

    catch
        eventTimes.Prem.wait_start = [];
    end
    eventTimes.Cor.cue = All{22, 1}{neuron, 1}.(eventCor);
    eventTimes.Cor.resp = All{22, 1}{neuron, 1}.(eventCorResp);
    eventTimes.Inc.cue = All{22, 1}{neuron, 1}.(eventInc);
    eventTimes.Inc.resp = All{22, 1}{neuron, 1}.(eventIncResp);
    eventTimes.Om.cue = All{22, 1}{neuron, 1}.(eventOm);
    eventTimes.Om.resp = All{22, 1}{neuron, 1}.(eventOm);
    try
        eventTimes.Prem.resp = All{22, 1}{neuron, 1}.(eventPremResp);
    catch
        eventTimes.Prem.resp = [];
    end
end

% Get binned values per trial
psthBinsValue.Cor.trial_start = [];
psthBinsValue.Inc.trial_start = [];
psthBinsValue.Om.trial_start = [];
psthBinsValue.Prem.trial_start = [];
psthBinsValue.Cor.wait_start = [];
psthBinsValue.Inc.wait_start = [];
psthBinsValue.Om.wait_start = [];
psthBinsValue.Prem.wait_start = [];
psthBinsValue.Cor.cue = [];
psthBinsValue.Inc.cue = [];
psthBinsValue.Om.cue = [];
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
try
    TrialDistr.PremTresholdExcluded = find(isnan(eventVideoBehavPrem.treshold));
    TrialDistr.PremTresholdIncluded = find(~isnan(eventVideoBehavPrem.treshold));
catch
    TrialDistr.PremTresholdExcluded = nan;
    TrialDistr.PremTresholdIncluded = nan;
end
    TrialDistr.CorTreshold = eventVideoBehavCor.treshold(TrialDistr.CorTresholdIncluded);
TrialDistr.IncTreshold = eventVideoBehavInc.treshold(TrialDistr.IncTresholdIncluded);
TrialDistr.OmTreshold = eventVideoBehavOm.treshold(TrialDistr.OmTresholdIncluded);
try
TrialDistr.PremTreshold = eventVideoBehavPrem.treshold(TrialDistr.PremTresholdIncluded);
catch
TrialDistr.PremTreshold = nan;
end

psthBinsValue.Cor.trial_start = binData(eventTimes.Cor.trial_start(TrialDistr.CorTresholdIncluded), All{21, 1}{neuron, 1}, TrialStart_tBef, TrialStart_tAft, binSize);
psthBinsValue.Inc.trial_start = binData(eventTimes.Inc.trial_start(TrialDistr.IncTresholdIncluded), All{21, 1}{neuron, 1}, TrialStart_tBef, TrialStart_tAft, binSize);
psthBinsValue.Om.trial_start = binData(eventTimes.Om.trial_start(TrialDistr.OmTresholdIncluded), All{21, 1}{neuron, 1}, TrialStart_tBef, TrialStart_tAft, binSize);
try
psthBinsValue.Prem.trial_start = binData(eventTimes.Prem.trial_start(TrialDistr.PremTresholdIncluded), All{21, 1}{neuron, 1}, TrialStart_tBef, TrialStart_tAft, binSize);
psthBinsValue.Prem.wait_start = binData(eventTimes.Prem.wait_start(TrialDistr.PremTresholdIncluded), All{21, 1}{neuron, 1}, Treshold_tBef, Treshold_tAft, binSize);
psthBinsValue.Prem.resp = binData(eventTimes.Prem.resp(TrialDistr.PremTresholdIncluded), All{21, 1}{neuron, 1}, Response_tBef, Response_tAft, binSize);
catch
    psthBinsValue.Prem.trial_start = nan;
psthBinsValue.Prem.wait_start = nan;
psthBinsValue.Prem.resp = nan;
end
psthBinsValue.Cor.wait_start = binData(eventTimes.Cor.wait_start(TrialDistr.CorTresholdIncluded), All{21, 1}{neuron, 1}, Treshold_tBef, Treshold_tAft, binSize);
psthBinsValue.Inc.wait_start = binData(eventTimes.Inc.wait_start(TrialDistr.IncTresholdIncluded), All{21, 1}{neuron, 1}, Treshold_tBef, Treshold_tAft, binSize);
psthBinsValue.Om.wait_start = binData(eventTimes.Om.wait_start(TrialDistr.OmTresholdIncluded), All{21, 1}{neuron, 1}, Treshold_tBef, Treshold_tAft, binSize);
psthBinsValue.Cor.resp = binData(eventTimes.Cor.resp(TrialDistr.CorTresholdIncluded), All{21, 1}{neuron, 1}, Response_tBef, Response_tAft, binSize);
psthBinsValue.Inc.resp = binData(eventTimes.Inc.resp(TrialDistr.IncTresholdIncluded), All{21, 1}{neuron, 1}, Response_tBef, Response_tAft, binSize);
psthBinsValue.Om.resp = binData(eventTimes.Om.resp(TrialDistr.OmTresholdIncluded), All{21, 1}{neuron, 1}, Response_tBef, Response_tAft, binSize);
psthBinsValue.Cor.cue = binData(eventTimes.Cor.cue(TrialDistr.CorTresholdIncluded), All{21, 1}{neuron, 1}, Cue_tBef, Cue_tAft, binSize);
psthBinsValue.Inc.cue = binData(eventTimes.Inc.cue(TrialDistr.IncTresholdIncluded), All{21, 1}{neuron, 1}, Cue_tBef, Cue_tAft, binSize);
psthBinsValue.Om.cue = binData(eventTimes.Om.cue(TrialDistr.OmTresholdIncluded), All{21, 1}{neuron, 1}, Cue_tBef, Cue_tAft, binSize);

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

nans = numel(find(isnan(psthBinsValue.Cor.cue(:,1))));
PSTH.Cor.cue = (nansum(psthBinsValue.Cor.cue,1)/(size(psthBinsValue.Cor.cue,1)-nans))*(1/binSize);
nans = numel(find(isnan(psthBinsValue.Inc.cue(:,1))));
PSTH.Inc.cue = (nansum(psthBinsValue.Inc.cue,1)/(size(psthBinsValue.Inc.cue,1)-nans))*(1/binSize);
nans = numel(find(isnan(psthBinsValue.Om.cue(:,1))));
PSTH.Om.cue = (nansum(psthBinsValue.Om.cue,1)/(size(psthBinsValue.Om.cue,1)-nans))*(1/binSize);

%% Make PSTH from randomly shuffled event times and calculate peak and information distance
switch InfoDist
    case 'yes'
        %Make randomized data
        outcome = {'Cor', 'Prem'};
        for n = 1:numel(outcome)

            timePeriod = {'TrialStart', 'Treshold','Cue','Response'};
            if strcmp(outcome{n}, 'Prem')
                timePeriod = {'TrialStart', 'Treshold', 'Response'};
            end
            for m = 1:numel(timePeriod)
                if strcmp(event,'ITI50') || strcmp(event,'All')
                    if strcmp(timePeriod{m},'TrialStart')
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 0;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 5;
                        eventType = 'trial_start';
                        tBef = 2;
                        tAft = 7;
                    elseif strcmp(timePeriod{m},'Treshold')
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 1;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 3;
                        eventType = 'wait_start';
                        tBef = 2;
                        tAft = 5;
                    elseif strcmp(timePeriod{m},'Cue') && ~strcmp(outcome{n}, 'Prem')
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 2;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 0;
                        eventType = 'cue';
                        tBef = 3;
                        tAft = 3;
                    elseif strcmp(timePeriod{m},'Response')
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 0;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 3;
                        eventType = 'resp';
                        tBef = 2;
                        tAft = 4;
                    end
                elseif strcmp(event,'ITI75')
                    if strcmp(timePeriod{m},'TrialStart')
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 0;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 7.5;
                        eventType = 'trial_start';
                        tBef = 2;
                        tAft = 10;
                    elseif strcmp(timePeriod{m},'Treshold')
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 1;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 3;
                        eventType = 'wait_start';
                        tBef = 2;
                        tAft = 8;
                    elseif strcmp(timePeriod{m},'Cue') && ~strcmp(outcome{n}, 'Prem')
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 2;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 0;
                        eventType = 'cue';
                        tBef = 2;
                        tAft = 2;
                    elseif strcmp(timePeriod{m},'Response')
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 0;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 3;
                        eventType = 'resp';
                        tBef = 2;
                        tAft = 4;
                    end
                elseif strcmp(event,'ITI125')
                    if strcmp(timePeriod{m},'TrialStart')
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 0;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 12.5;
                        eventType = 'trial_start';
                        tBef = 2;
                        tAft = 15;
                    elseif strcmp(timePeriod{m},'Treshold')
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 1;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 3;
                        eventType = 'wait_start';
                        tBef = 2;
                        tAft = 13;
                    elseif strcmp(timePeriod{m},'Cue') && ~strcmp(outcome{n}, 'Prem')
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 2;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 0;
                        eventType = 'cue';
                        tBef = 3;
                        tAft = 3;
                    elseif strcmp(timePeriod{m},'Response')
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID = 0;
                        PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID = 3;
                        eventType = 'resp';
                        tBef = 2;
                        tAft = 4;
                    end
                end
            
                binSizeID = 1;
            if ~(strcmp(timePeriod{m},'Cue') && strcmp(outcome{n}, 'Prem'))
            BinsPerSec = 1/binSizeID;
            FRmean = mean(PSTH.Cor.trial_start((TrialStart_tBef*(1/binSizeID)-1*(1/binSizeID)):TrialStart_tBef*(1/binSizeID)));
            PSTH.ID.(outcome{n}).(timePeriod{m}).ID_real = InformDistance(FRmean, PSTH.(outcome{n}).(eventType), tBef, PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID, PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID, binSizeID);
            %maxPeak = max(smoothdata(PSTH.(outcome{n}).(eventType)((tBef*BinsPerSec)-(PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID*BinsPerSec):(tBef*BinsPerSec)+(PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID*BinsPerSec))-All{15,1}(neuron,2),'gaussian',3));
            %minPeak = min(smoothdata(PSTH.(outcome{n}).(eventType)((tBef*BinsPerSec)-(PSTH.ID.(outcome{n}).(timePeriod{m}).tBef_ID*BinsPerSec):(tBef*BinsPerSec)+(PSTH.ID.(outcome{n}).(timePeriod{m}).tAft_ID*BinsPerSec))-All{15,1}(neuron,2),'gaussian',3));
            PSTH.ID.(outcome{n}).(timePeriod{m}).Peak_real = max(abs([minPeak maxPeak]));
            if PSTH.ID.(outcome{n}).(timePeriod{m}).Peak_real == abs(minPeak)
                PSTH.ID.(outcome{n}).(timePeriod{m}).Peak_PosOrNeg = -1;
            elseif PSTH.ID.(outcome{n}).(timePeriod{m}).Peak_real == abs(maxPeak)
                PSTH.ID.(outcome{n}).(timePeriod{m}).Peak_PosOrNeg = 1;
            end

            end
            end
        end
    case 'no'
end


%% Plotting

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

            lightIxs = find(All{1,1}{neuron,1}(:,3) < 0.01 & All{1,1}{neuron,1}(:,15)>= parameters_OI.PWCor & All{1,1}{neuron,1}(:,2)>= parameters_OI.minLightInt & All{1,1}{neuron,1}(:,6)>= parameters_OI.ReliabilityTreshold);                %stimBlocksPower = find(All{1,1}{neuron,1}(lightIxs,2) == 100 | All{1,1}{neuron,1}(lightIxs,2) == 75);
            [~, stimBlocksDur] = max(All{1,1}{neuron,1}(lightIxs,1));
            stimBlock = lightIxs(stimBlocksDur);

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

            figure

                dim = [0.7 0.15 0.3 0.3];
                content = {sprintf('Avg Spikes s^{-1} = %g',All{15,1}(neuron,2)), sprintf('ID %g', All{3, 1}(neuron,1)), sprintf('L-ratio %g', All{3, 1}(neuron,2)), sprintf('ISI1.5 %g',...
                All{3, 1}(neuron,7)), sprintf('session dur (hour) %g', All{22, 1}{neuron, 1}.corResp(end)/3600), sprintf('shank %g', All{4, 1}(neuron,1)),...
                sprintf('depth %g', All{11, 1}(neuron,1)), sprintf('animal %s', All{8, 1}{neuron,1}), sprintf('buddy %g', All{25, 1}(neuron,1)), sprintf('SALT_P %g', min(All{1, 1}{neuron,1}(:,3)))};
                annotation('textbox',dim,'String',content,'FitBoxToText','on');
                plotedit on
            
            % Get timeaxis for PSTH
            time_axis_Treshold = round(linspace(-Treshold_tBef+(binSize),Treshold_tAft,((Treshold_tBef+Treshold_tAft)/binSize)),1);
            time_axis_TrialStart = round(linspace(-TrialStart_tBef+(binSize),TrialStart_tAft,((TrialStart_tBef+TrialStart_tAft)/binSize)),1);
            time_axis_Response = round(linspace(-Response_tBef+(binSize),Response_tAft,((Response_tBef+Response_tAft)/binSize)),1);
            time_axis_Cue = round(linspace(-Cue_tBef+(binSize),Cue_tAft,((Cue_tBef+Cue_tAft)/binSize)),1);

            subplot(3,4,5)
            hold on
            SEM = nanstd(psthBinsValue.Cor.trial_start,1)/sqrt(size(psthBinsValue.Cor.trial_start,1)); 
            shadedErrorBar(time_axis_TrialStart,smooth(PSTH.Cor.trial_start,3),SEM,{'k'});
            SEM = nanstd(psthBinsValue.Prem.trial_start,1)/sqrt(size(psthBinsValue.Prem.trial_start,1)); 
            shadedErrorBar(time_axis_TrialStart,smooth(PSTH.Prem.trial_start,3),SEM,{'c'});
            line([-TrialStart_tBef TrialStart_tAft],[mean(PSTH.Cor.trial_start(1:10)) mean(PSTH.Cor.trial_start(1:10))])

            axis([-TrialStart_tBef TrialStart_tAft 0 inf])

            subplot(3,4,6)
            hold on
            SEM = nanstd(psthBinsValue.Cor.wait_start,1)/sqrt(size(psthBinsValue.Cor.wait_start,1)); 
            shadedErrorBar(time_axis_Treshold,smooth(PSTH.Cor.wait_start,3),SEM,{'k'});
            SEM = nanstd(psthBinsValue.Prem.wait_start,1)/sqrt(size(psthBinsValue.Prem.wait_start,1)); 
            shadedErrorBar(time_axis_Treshold,smooth(PSTH.Prem.wait_start,3),SEM,{'c'});
            axis([-Treshold_tBef Treshold_tAft 0 inf])

            subplot(3,4,7)
            hold on
            SEM = nanstd(psthBinsValue.Cor.cue,1)/sqrt(size(psthBinsValue.Cor.cue,1)); 
            shadedErrorBar(time_axis_Cue,smooth(PSTH.Cor.cue,3),SEM,{'k'});
            SEM = nanstd(psthBinsValue.Prem.resp(:,1:20),1)/sqrt(size(psthBinsValue.Prem.resp(:,1:20),1)); 
            shadedErrorBar(time_axis_Cue,smooth(PSTH.Prem.resp(:,1:20),3),SEM,{'c'});
            axis([-Cue_tBef Cue_tAft 0 inf])

            subplot(3,4,8)
            hold on
            SEM = nanstd(psthBinsValue.Cor.resp,1)/sqrt(size(psthBinsValue.Cor.resp,1)); 
            shadedErrorBar(time_axis_Response,smooth(PSTH.Cor.resp,3),SEM,{'k'});
            SEM = nanstd(psthBinsValue.Prem.resp,1)/sqrt(size(psthBinsValue.Prem.resp,1)); 
            shadedErrorBar(time_axis_Response,smooth(PSTH.Prem.resp,3),SEM,{'c'});
            axis([-Response_tBef Response_tAft 0 inf])
            
            dotSize = 0.1;
            subplot(3,4,9);
            hold on
            n = 1;
            for i=1:numel(TrialDistr.CorTresholdIncluded)
                spnum=find(All{21, 1}{neuron, 1}>eventTimes.Cor.trial_start(TrialDistr.CorTresholdIncluded(i))-TrialStart_tBef & All{21, 1}{neuron, 1} < eventTimes.Cor.trial_start(TrialDistr.CorTresholdIncluded(i))+TrialStart_tAft);
                scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Cor.trial_start(TrialDistr.CorTresholdIncluded(i)),repmat(n,numel(spnum),1)',dotSize,'k','filled');
                scatter(eventVideoBehavCor.treshold(i),repmat(n,numel(1),1)',dotSize*2,'m','filled');
                scatter(eventTimes.Cor.resp(TrialDistr.CorTresholdIncluded(i))-eventTimes.Cor.trial_start(TrialDistr.CorTresholdIncluded(i)),repmat(n,numel(1),1)',dotSize*2,'c','filled');

                n = n+1;
            end
            for i=1:numel(TrialDistr.PremTresholdIncluded)
                spnum=find(All{21, 1}{neuron, 1}>eventTimes.Prem.trial_start(TrialDistr.PremTresholdIncluded(i))-TrialStart_tBef & All{21, 1}{neuron, 1} < eventTimes.Prem.trial_start(TrialDistr.PremTresholdIncluded(i))+TrialStart_tAft);
                scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Prem.trial_start(TrialDistr.PremTresholdIncluded(i)),repmat(n,numel(spnum),1)',dotSize,[0.8500 0.3250 0.0980],'filled');
                scatter(eventVideoBehavPrem.treshold(i),repmat(n,numel(1),1)',dotSize*2,'m','filled');
                scatter(eventTimes.Prem.resp(TrialDistr.PremTresholdIncluded(i))-eventTimes.Prem.trial_start(TrialDistr.PremTresholdIncluded(i)),repmat(n,numel(1),1)',dotSize*2,'c','filled');

                n = n+1;
            end
            axis([-TrialStart_tBef TrialStart_tAft 1 numel(TrialDistr.CorTresholdIncluded)+numel(TrialDistr.PremTresholdIncluded)])
            ylabel('Trial')
            xlabel('time (sec)')

            subplot(3,4,10);
            hold on
            n = 1;
            for i=1:numel(TrialDistr.CorTresholdIncluded)
                spnum=find(All{21, 1}{neuron, 1}>eventTimes.Cor.wait_start(TrialDistr.CorTresholdIncluded(i))-Treshold_tBef & All{21, 1}{neuron, 1} < eventTimes.Cor.wait_start(TrialDistr.CorTresholdIncluded(i))+Treshold_tAft);
                scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Cor.wait_start(TrialDistr.CorTresholdIncluded(i)),repmat(n,numel(spnum),1)',dotSize,'k','filled');
                n = n+1;
            end
            for i=1:numel(TrialDistr.PremTresholdIncluded)
                spnum=find(All{21, 1}{neuron, 1}>eventTimes.Prem.wait_start(TrialDistr.PremTresholdIncluded(i))-Treshold_tBef & All{21, 1}{neuron, 1} < eventTimes.Prem.wait_start(TrialDistr.PremTresholdIncluded(i))+Treshold_tAft);
                scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Prem.wait_start(TrialDistr.PremTresholdIncluded(i)),repmat(n,numel(spnum),1)',dotSize,[0.8500 0.3250 0.0980],'filled');
                n = n+1;
            end
            axis([-Treshold_tBef Treshold_tAft 1 numel(TrialDistr.CorTresholdIncluded)+numel(TrialDistr.PremTresholdIncluded)])
            ylabel('Trial')
            xlabel('time (sec)')

            subplot(3,4,11);
            hold on
            n = 1;
            for i=1:numel(eventTimes.Cor.cue)
                spnum=find(All{21, 1}{neuron, 1}>eventTimes.Cor.cue(i)-Response_tBef & All{21, 1}{neuron, 1} < eventTimes.Cor.cue(i)+Response_tAft);
                scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Cor.cue(i),repmat(n,numel(spnum),1)',dotSize,'k','filled');
                n = n+1;
            end
            for i=1:numel(eventTimes.Prem.resp)
                spnum=find(All{21, 1}{neuron, 1}>eventTimes.Prem.resp(i)-Response_tBef & All{21, 1}{neuron, 1} < eventTimes.Prem.resp(i)+Response_tAft);
                scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Prem.resp(i),repmat(n,numel(spnum),1)',dotSize,[0.8500 0.3250 0.0980],'filled');
                n = n+1;
            end
            axis([-Response_tBef Response_tAft 1 numel([eventTimes.Cor.cue eventTimes.Prem.resp])])
            ylabel('Trial')
            xlabel('time (sec)')

            subplot(3,4,12);
            hold on
            n = 1;
            for i=1:numel(TrialDistr.CorTresholdIncluded)
                spnum=find(All{21, 1}{neuron, 1}>eventTimes.Cor.resp(TrialDistr.CorTresholdIncluded(i))-Response_tBef & All{21, 1}{neuron, 1} < eventTimes.Cor.resp(TrialDistr.CorTresholdIncluded(i))+Response_tAft);
                scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Cor.resp(TrialDistr.CorTresholdIncluded(i)),repmat(n,numel(spnum),1)',dotSize,'k','filled');
                n = n+1;
            end
            for i=1:numel(TrialDistr.PremTresholdIncluded)
                spnum=find(All{21, 1}{neuron, 1}>eventTimes.Prem.resp(TrialDistr.PremTresholdIncluded(i))-Response_tBef & All{21, 1}{neuron, 1} < eventTimes.Prem.resp(TrialDistr.PremTresholdIncluded(i))+Response_tAft);
                scatter(All{21, 1}{neuron, 1}(spnum)-eventTimes.Prem.resp(TrialDistr.PremTresholdIncluded(i)),repmat(n,numel(spnum),1)',dotSize,[0.8500 0.3250 0.0980],'filled');
                n = n+1;
            end
            axis([-Response_tBef Response_tAft 1 numel(TrialDistr.CorTresholdIncluded)+numel(TrialDistr.PremTresholdIncluded)])
            ylabel('Trial')
            xlabel('time (sec)')

    end
end
end
