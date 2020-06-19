function [VideoTrialsRemoval, eventVideoBehavCor, eventVideoBehavInc, eventVideoBehavOm] = videoTrialRemoval_SD(All, txtFile, event)
 

    Ixes = find(contains(All{5,1}(:,2),txtFile));
    neuron = Ixes(1);


    if strcmp(event,'All')
        eventCor = 'corCue';
        eventInc = 'incCue';
        eventOm = 'omCue';
    elseif strcmp(event,'SD2')
        eventCor = 'corSD2Cue';
        eventInc = 'incSD2Cue';
        eventOm = 'omSD2Cue';
    elseif strcmp(event,'SD5')
        eventCor = 'corSD5Cue';
        eventInc = 'incSD5Cue';
        eventOm = 'omSD5Cue';
    elseif strcmp(event,'SD10')
        eventCor = 'corSD10Cue';
        eventInc = 'incSD10Cue';
        eventOm = 'omSD10Cue';
    end
    
    
    eventVideoBehavCor = eventVideo(All, eventCor, neuron, 0, 5,'PSTH_EVENT');
    eventVideoBehavInc = eventVideo(All, eventInc, neuron, 0, 5,'PSTH_EVENT');
    eventVideoBehavOm = eventVideo(All, eventOm, neuron, 0, 5,'PSTH_EVENT');
    
    VideoTrialsRemoval.Cor_ThreeSecsIx = find(eventVideoBehavCor.treshold <= 3 & eventVideoBehavCor.tresholdSecondPass == 0 & eventVideoBehavCor.tresholdLateralDist >= 0.03);
    VideoTrialsRemoval.Inc_ThreeSecsIx = find(eventVideoBehavInc.treshold <= 3 & eventVideoBehavInc.tresholdSecondPass == 0 & eventVideoBehavInc.tresholdLateralDist >= 0.03);
    VideoTrialsRemoval.Om_ThreeSecsIx = find(eventVideoBehavOm.treshold <= 3 & eventVideoBehavOm.tresholdSecondPass == 0 & eventVideoBehavOm.tresholdLateralDist >= 0.03);
    VideoTrialsRemoval.Cor_ThreeSecsIxRemoved = find(eventVideoBehavCor.treshold > 3 | eventVideoBehavCor.tresholdSecondPass == 1 | eventVideoBehavCor.tresholdLateralDist < 0.03);
    VideoTrialsRemoval.Inc_ThreeSecsIxRemoved = find(eventVideoBehavInc.treshold > 3 | eventVideoBehavInc.tresholdSecondPass == 1 | eventVideoBehavInc.tresholdLateralDist < 0.03);
    VideoTrialsRemoval.Om_ThreeSecsIxRemoved = find(eventVideoBehavOm.treshold > 3 | eventVideoBehavOm.tresholdSecondPass == 1 | eventVideoBehavOm.tresholdLateralDist < 0.03);
    