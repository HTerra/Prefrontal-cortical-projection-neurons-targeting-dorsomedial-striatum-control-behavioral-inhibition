function All = insertBehavVars(All)

for NeuPos = 1:size(All{22, 1},1)
    NeuPos
    CueOn = [All{22, 1}{NeuPos,1}.CH1 All{22, 1}{NeuPos,1}.CH2 All{22, 1}{NeuPos,1}.CH3 All{22, 1}{NeuPos,1}.CH4 All{22, 1}{NeuPos,1}.CH5];
    for i = 1:size(All{22, 1}{NeuPos,1}.CH6,2)
        [cCor(i), indexCor(i)] = min(abs(CueOn-All{22, 1}{NeuPos,1}.CH6(i)));
    end
    for i = 1:size(All{22, 1}{NeuPos,1}.CH8,2)
        [cInc(i), indexInc(i)] = min(abs(CueOn-All{22, 1}{NeuPos,1}.CH8(i)));
    end
    for i = 1:size(All{22, 1}{NeuPos,1}.CH10,2)
        [cOm(i), indexOm(i)] = min(abs(CueOn-All{22, 1}{NeuPos,1}.CH10(i)));
    end
    for i = 1:size(All{22, 1}{NeuPos,1}.CH12,2)
        [cPrem(i), indexPrem(i)] = min(abs(CueOn-All{22, 1}{NeuPos,1}.CH12(i)));
    end
    for i = 1:size(All{22, 1}{NeuPos,1}.CH13,2)
        [cPers(i), indexPers(i)] = min(abs(CueOn-All{22, 1}{NeuPos,1}.CH13(i)));
    end
    
    All{22, 1}{NeuPos,1}.CH17 = All{22, 1}{NeuPos,1}.CH6-cCor;
    All{22, 1}{NeuPos,1}.CH18 = All{22, 1}{NeuPos,1}.CH8-cInc;
    All{22, 1}{NeuPos,1}.CH19 = All{22, 1}{NeuPos,1}.CH10-cOm;
    All{22, 1}{NeuPos,1}.CH20 = All{22, 1}{NeuPos,1}.CH12-cPrem;
    All{22, 1}{NeuPos,1}.CH21 = All{22, 1}{NeuPos,1}.CH13-cPers;
    
    clear cCor indexCor cInc indexInc cOm indexOm cPrem indexPrem cPers indexPers
    
    % Rename struct fields to logical names
    All{22, 1}{NeuPos,1}.corResp = All{22, 1}{NeuPos,1}.CH6;
    All{22, 1}{NeuPos,1}.incResp = All{22, 1}{NeuPos,1}.CH8;
    All{22, 1}{NeuPos,1}.omResp = All{22, 1}{NeuPos,1}.CH10;
    All{22, 1}{NeuPos,1}.premResp = All{22, 1}{NeuPos,1}.CH12;
    All{22, 1}{NeuPos,1}.persResp = All{22, 1}{NeuPos,1}.CH13;
    
    All{22, 1}{NeuPos,1}.corCue = All{22, 1}{NeuPos,1}.CH17;
    All{22, 1}{NeuPos,1}.incCue = All{22, 1}{NeuPos,1}.CH18;
    All{22, 1}{NeuPos,1}.omCue = All{22, 1}{NeuPos,1}.CH19;
    
    if strcmp(All{5, 1} {NeuPos,1},'VAR_ITI')
        [Correct, Incorrect, Omission, Premature, Perserverative, nrTrialsStage, perf{NeuPos,1}]  = CC_ephys_var_ITI_SD_TTL_Master(10,1,All{5, 1} {NeuPos,2});
        
        % TEMP FIX LOOK AK WHY THERE IS MISMATCH IN AMOUNT OF CORRECT AND
        % PREMATURE TRIALS. PROBABLY FEW MISMATCH TRIAL TYPES
%         if ismember(NeuPos,232:246)
%             Correct.TrialPosITI50(74) = [];
%             Premature.TrialPosITI125(22:23) = [];
%         end
        
%         if ismember(NeuPos,176:190)
%             Correct.TrialPosITI50(74) = [];
%             Premature.TrialPosITI125(22:23) = [];
%             Correct.MagLatency = Correct.MagLatency(1:end-1);
%             All{22, 1}{NeuPos,1}.corMag = Correct.MagLatency + All{22, 1}{NeuPos,1}.corResp;
%             Correct.MagLatencyITI50 = Correct.MagLatencyITI50(1:end-1);
%         end
        
        % Get premature and corMagLatency and TTL time
        All{22, 1}{NeuPos,1}.corMagLatency = Correct.MagLatency;
        All{22, 1}{NeuPos,1}.corMag = All{22, 1}{NeuPos,1}.corMagLatency+All{22, 1}{NeuPos,1}.corResp;
        All{22, 1}{NeuPos,1}.preLatency = Premature.Latency;
        
        % Get timestamps for VAR ITI and VAR SD conditions
        All{22, 1}{NeuPos,1}.corITI50Resp = All{22, 1}{NeuPos,1}.CH6(Correct.TrialPosITI50);
        All{22, 1}{NeuPos,1}.corITI75Resp = All{22, 1}{NeuPos,1}.CH6(Correct.TrialPosITI75);
        All{22, 1}{NeuPos,1}.corITI125Resp = All{22, 1}{NeuPos,1}.CH6(Correct.TrialPosITI125);
        All{22, 1}{NeuPos,1}.incITI50Resp = All{22, 1}{NeuPos,1}.CH8(Incorrect.TrialPosITI50);
        All{22, 1}{NeuPos,1}.incITI75Resp = All{22, 1}{NeuPos,1}.CH8(Incorrect.TrialPosITI75);
        All{22, 1}{NeuPos,1}.incITI125Resp = All{22, 1}{NeuPos,1}.CH8(Incorrect.TrialPosITI125);
        All{22, 1}{NeuPos,1}.omITI50Resp = All{22, 1}{NeuPos,1}.CH10(Omission.TrialPosITI50);
        All{22, 1}{NeuPos,1}.omITI75Resp = All{22, 1}{NeuPos,1}.CH10(Omission.TrialPosITI75);
        All{22, 1}{NeuPos,1}.omITI125Resp = All{22, 1}{NeuPos,1}.CH10(Omission.TrialPosITI125);
        All{22, 1}{NeuPos,1}.premITI50Resp = All{22, 1}{NeuPos,1}.CH12(Premature.TrialPosITI50);
        All{22, 1}{NeuPos,1}.premITI75Resp = All{22, 1}{NeuPos,1}.CH12(Premature.TrialPosITI75);
        All{22, 1}{NeuPos,1}.premITI125Resp = All{22, 1}{NeuPos,1}.CH12(Premature.TrialPosITI125);
        All{22, 1}{NeuPos,1}.persITI50Resp = All{22, 1}{NeuPos,1}.CH13(Perserverative.TrialPosITI50);
        All{22, 1}{NeuPos,1}.persITI75Resp = All{22, 1}{NeuPos,1}.CH13(Perserverative.TrialPosITI75);
        All{22, 1}{NeuPos,1}.persITI125Resp = All{22, 1}{NeuPos,1}.CH13(Perserverative.TrialPosITI125);

        All{22, 1}{NeuPos,1}.corITI50Cue = All{22, 1}{NeuPos,1}.CH17(Correct.TrialPosITI50);
        All{22, 1}{NeuPos,1}.corITI75Cue = All{22, 1}{NeuPos,1}.CH17(Correct.TrialPosITI75);
        All{22, 1}{NeuPos,1}.corITI125Cue = All{22, 1}{NeuPos,1}.CH17(Correct.TrialPosITI125);
        All{22, 1}{NeuPos,1}.incITI50Cue = All{22, 1}{NeuPos,1}.CH18(Incorrect.TrialPosITI50);
        All{22, 1}{NeuPos,1}.incITI75Cue = All{22, 1}{NeuPos,1}.CH18(Incorrect.TrialPosITI75);
        All{22, 1}{NeuPos,1}.incITI125Cue = All{22, 1}{NeuPos,1}.CH18(Incorrect.TrialPosITI125);
        All{22, 1}{NeuPos,1}.omITI50Cue = All{22, 1}{NeuPos,1}.CH19(Omission.TrialPosITI50);
        All{22, 1}{NeuPos,1}.omITI75Cue = All{22, 1}{NeuPos,1}.CH19(Omission.TrialPosITI75);
        All{22, 1}{NeuPos,1}.omITI125Cue = All{22, 1}{NeuPos,1}.CH19(Omission.TrialPosITI125);
       
        All{22, 1}{NeuPos,1}.corMagLatencyITI5 = Correct.MagLatencyITI50;
		All{22, 1}{NeuPos,1}.corMagLatencyITI75 = Correct.MagLatencyITI75;
		All{22, 1}{NeuPos,1}.corMagLatencyITI125 = Correct.MagLatencyITI125;
        All{22, 1}{NeuPos,1}.corMagITI5 = Correct.MagLatencyITI50 + All{22, 1}{NeuPos,1}.corITI50Resp;
		All{22, 1}{NeuPos,1}.corMagITI75 = Correct.MagLatencyITI75 + All{22, 1}{NeuPos,1}.corITI75Resp;
		All{22, 1}{NeuPos,1}.corMagITI125 = Correct.MagLatencyITI125 + All{22, 1}{NeuPos,1}.corITI125Resp;
        
        All{22, 1}{NeuPos,1}.preLatencyITI5 = Premature.LatencyITI50;
		All{22, 1}{NeuPos,1}.preLatencyITI75 = Premature.LatencyITI75;
		All{22, 1}{NeuPos,1}.preLatencyITI125 = Premature.LatencyITI125;
        
        All{22, 1}{NeuPos,1}.preTrialStartITI5 = All{22, 1}{NeuPos,1}.premITI50Resp-All{22, 1}{NeuPos,1}.preLatencyITI5;
		All{22, 1}{NeuPos,1}.preTrialStartITI75 = All{22, 1}{NeuPos,1}.premITI75Resp-All{22, 1}{NeuPos,1}.preLatencyITI75;
		All{22, 1}{NeuPos,1}.preTrialStartITI125 = All{22, 1}{NeuPos,1}.premITI125Resp-All{22, 1}{NeuPos,1}.preLatencyITI125;
        
    elseif strcmp(All{5, 1} {NeuPos,1},'VAR_SD')
        [Correct, Incorrect, Omission, Premature, Perserverative, nrTrialsStage, perf{NeuPos,1}]  = CC_ephys_var_ITI_SD_TTL_Master(10,1,All{5, 1} {NeuPos,2});
        
        
        % TEMP FIX LOOK AK WHY THERE IS MISMATCH IN AMOUNT OF CORRECT AND
        % PREMATURE TRIALS. PROBABLY FEW MISMATCH TRIAL TYPES
%         if ismember(NeuPos,323:340)
%             Correct.TrialPosSD10(143) = [];
%             Correct.MagLatency = Correct.MagLatency(1:end-1);
%             All{22, 1}{NeuPos,1}.corMag = Correct.MagLatency + All{22, 1}{NeuPos,1}.corResp;
%         end
        
        % Get premature and corMagLatency and TTL time
        All{22, 1}{NeuPos,1}.corMagLatency = Correct.MagLatency;
        All{22, 1}{NeuPos,1}.corMag = Correct.MagLatency+All{22, 1}{NeuPos,1}.corResp;
        All{22, 1}{NeuPos,1}.preLatency = Premature.Latency;
        
        All{22, 1}{NeuPos,1}.corSD2Resp = All{22, 1}{NeuPos,1}.CH6(Correct.TrialPosSD2);
        All{22, 1}{NeuPos,1}.corSD5Resp = All{22, 1}{NeuPos,1}.CH6(Correct.TrialPosSD5);
        All{22, 1}{NeuPos,1}.corSD10Resp = All{22, 1}{NeuPos,1}.CH6(Correct.TrialPosSD10);
        All{22, 1}{NeuPos,1}.incSD2Resp = All{22, 1}{NeuPos,1}.CH8(Incorrect.TrialPosSD2);
        All{22, 1}{NeuPos,1}.incSD5Resp = All{22, 1}{NeuPos,1}.CH8(Incorrect.TrialPosSD5);
        All{22, 1}{NeuPos,1}.incSD10Resp = All{22, 1}{NeuPos,1}.CH8(Incorrect.TrialPosSD10);
        All{22, 1}{NeuPos,1}.omSD2Resp = All{22, 1}{NeuPos,1}.CH10(Omission.TrialPosSD2);
        All{22, 1}{NeuPos,1}.omSD5Resp = All{22, 1}{NeuPos,1}.CH10(Omission.TrialPosSD5);
        All{22, 1}{NeuPos,1}.omSD10Resp = All{22, 1}{NeuPos,1}.CH10(Omission.TrialPosSD10);
        All{22, 1}{NeuPos,1}.premSD2Resp = All{22, 1}{NeuPos,1}.CH12(Premature.TrialPosSD2);
        All{22, 1}{NeuPos,1}.premSD5Resp = All{22, 1}{NeuPos,1}.CH12(Premature.TrialPosSD5);
        All{22, 1}{NeuPos,1}.premSD10Resp = All{22, 1}{NeuPos,1}.CH12(Premature.TrialPosSD10);
        All{22, 1}{NeuPos,1}.persSD2Resp = All{22, 1}{NeuPos,1}.CH13(Perserverative.TrialPosSD2);
        All{22, 1}{NeuPos,1}.persSD5Resp = All{22, 1}{NeuPos,1}.CH13(Perserverative.TrialPosSD5);
        All{22, 1}{NeuPos,1}.persSD10Resp = All{22, 1}{NeuPos,1}.CH13(Perserverative.TrialPosSD10);

        All{22, 1}{NeuPos,1}.corSD2Cue = All{22, 1}{NeuPos,1}.CH17(Correct.TrialPosSD2);
        All{22, 1}{NeuPos,1}.corSD5Cue = All{22, 1}{NeuPos,1}.CH17(Correct.TrialPosSD5);
        All{22, 1}{NeuPos,1}.corSD10Cue = All{22, 1}{NeuPos,1}.CH17(Correct.TrialPosSD10);
        All{22, 1}{NeuPos,1}.incSD2Cue = All{22, 1}{NeuPos,1}.CH18(Incorrect.TrialPosSD2);
        All{22, 1}{NeuPos,1}.incSD5Cue = All{22, 1}{NeuPos,1}.CH18(Incorrect.TrialPosSD5);
        All{22, 1}{NeuPos,1}.incSD10Cue = All{22, 1}{NeuPos,1}.CH18(Incorrect.TrialPosSD10);
        All{22, 1}{NeuPos,1}.omSD2Cue = All{22, 1}{NeuPos,1}.CH19(Omission.TrialPosSD2);
        All{22, 1}{NeuPos,1}.omSD5Cue = All{22, 1}{NeuPos,1}.CH19(Omission.TrialPosSD5);
        All{22, 1}{NeuPos,1}.omSD10Cue = All{22, 1}{NeuPos,1}.CH19(Omission.TrialPosSD10);
        
                
        All{22, 1}{NeuPos,1}.corMagLatencySD2 = Correct.MagLatencySD2;
		All{22, 1}{NeuPos,1}.corMagLatencySD5 = Correct.MagLatencySD5;
		All{22, 1}{NeuPos,1}.corMagLatencySD10 = Correct.MagLatencySD10;
        All{22, 1}{NeuPos,1}.corMagSD2 = Correct.MagLatencySD2 + All{22, 1}{NeuPos,1}.corSD2Resp;
		All{22, 1}{NeuPos,1}.corMagSD5 = Correct.MagLatencySD5 + All{22, 1}{NeuPos,1}.corSD5Resp;
		All{22, 1}{NeuPos,1}.corMagSD10 = Correct.MagLatencySD10 + All{22, 1}{NeuPos,1}.corSD10Resp;
        
        All{22, 1}{NeuPos,1}.preLatencySD2 = Premature.LatencySD2;
		All{22, 1}{NeuPos,1}.preLatencySD5 = Premature.LatencySD5;
		All{22, 1}{NeuPos,1}.preLatencySD10 = Premature.LatencySD10;
    else
        [Correct, Incorrect, Omission, Premature, Perserverative, nrTrialsStage, perf{NeuPos,1}]  = CC_ephys_var_ITI_SD_TTL_Master(9,1,All{5, 1} {NeuPos,2});
    end
    
    %TimeCues = [All{22, 1}{NeuPos,1}.CH1 All{22, 1}{NeuPos,1}.CH2 All{22, 1}{NeuPos,1}.CH3 All{22, 1}{NeuPos,1}.CH4 All{22, 1}{NeuPos,1}.CH5];
    %TimeLastCue = TimeCues(end);
    %All{22, 1}{NeuPos,1}.TrialStartEnd(NeuPos,1) = 0;
    %All{22, 1}{NeuPos,1}.TrialStartEnd(NeuPos,2) = TimeCues(end);
    
    %SpikeIxRem = All{21, 1}{NeuPos,1}>(TimeLastCue+30);
    %All{21, 1}{NeuPos,1}(find(SpikeIxRem)) = [];
    
    
    %All{22, 1}{NeuPos,1} = rmfield(All{22, 1}{NeuPos,1},{'CH6', 'CH8', 'CH10', 'CH12' ,'CH13' ,'CH17','CH18','CH19'});
end