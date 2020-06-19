function [VideoTrialsRemoval, eventVideoBehavCor, eventVideoBehavPrem] = videoTrialRemoval_ITI(All, txtFile, event)
 

    Ixes = find(contains(All{5,1}(:,2),txtFile));
    neuron = Ixes(1);


    if strcmp(event,'All')
        eventCor = 'corCue';
        eventCorResp = 'corResp';
        eventPremResp = 'premResp';
    elseif strcmp(event,'ITI50')
        eventCor = 'corITI50Cue';
        eventPremResp = 'premITI50Resp';
        eventTimes.Cor.trial_start = All{22, 1}{neuron, 1}.(eventCor)-5;
        eventTimes.Prem.trial_start = All{22, 1}{neuron, 1}.(eventPremResp)-All{22, 1}{neuron, 1}.preLatencyITI5;
        eventVideoBehavCor = eventVideo(All, eventCor, neuron, 0, 5,'PSTH_EVENT');
        eventVideoBehavPrem = eventVideo(All, eventPremResp, neuron, 0, 5,'PSTH_EVENT');
    elseif strcmp(event,'ITI75')
        eventCor = 'corITI75Cue';
        eventPremResp = 'premITI75Resp';
        eventTimes.Cor.trial_start = All{22, 1}{neuron, 1}.(eventCor)-7.5;
        eventTimes.Prem.trial_start = All{22, 1}{neuron, 1}.(eventPremResp)-All{22, 1}{neuron, 1}.preLatencyITI75;
        eventVideoBehavCor = eventVideo(All, eventCor, neuron, 0, 7.5,'PSTH_EVENT');
        eventVideoBehavPrem = eventVideo(All, eventPremResp, neuron, 0, 7.5,'PSTH_EVENT');
    elseif strcmp(event,'ITI125')
        eventCor = 'corITI125Cue';
        eventPremResp = 'premITI125Resp';
        eventTimes.Cor.trial_start = All{22, 1}{neuron, 1}.(eventCor)-12.5;
        eventTimes.Prem.trial_start = All{22, 1}{neuron, 1}.(eventPremResp)-All{22, 1}{neuron, 1}.preLatencyITI125;
        eventVideoBehavCor = eventVideo(All, eventCor, neuron, 0, 12.5,'PSTH_EVENT');
        eventVideoBehavPrem = eventVideo(All, eventPremResp, neuron, 0, 12.5,'PSTH_EVENT');
    end
    
    
    
    %% Remove trials where animal went over treshold later than three seconds, passed the treshold again after first pass or had a lateral distance covered lower than 0.03 after treshold crossing
    
    try
        %VideoTrialsRemoval.Prem_ThreeSecsIx = find(eventVideoBehavPrem.treshold <= 3 & eventVideoBehavPrem.tresholdSecondPass == 0 & eventVideoBehavPrem.tresholdLateralDist >= 0.03);
        VideoTrialsRemoval.Prem_ThreeSecsIx = find(eventVideoBehavPrem.tresholdSecondPass == 0);
    catch
        VideoTrialsRemoval.Prem_ThreeSecsIx = [];
    end
    try
        %VideoTrialsRemoval.Prem_ThreeSecsIxRemoved = find(eventVideoBehavPrem.treshold > 3 | eventVideoBehavPrem.tresholdSecondPass == 1 | eventVideoBehavPrem.tresholdLateralDist < 0.03);
        VideoTrialsRemoval.Prem_ThreeSecsIxRemoved = find(eventVideoBehavPrem.tresholdSecondPass == 1);

    catch
        VideoTrialsRemoval.Prem_ThreeSecsIxRemoved = [];
    end