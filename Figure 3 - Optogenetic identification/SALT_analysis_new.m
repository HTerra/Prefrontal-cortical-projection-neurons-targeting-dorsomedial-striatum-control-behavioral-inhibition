for neuron = 1:size(All{1, 1},1)
    for stimBlock = 1:size(All{1, 1}{neuron, 1},1)
        [spt_test, spt_baseline, FSLatency, jitter, reliability, spt_spikeIx, baselineTimeArray] = binMakerSALT2_FORSALTREANALYSIS(All, neuron, stimBlock);
        [p I] = SALT(spt_baseline,spt_test,0.001,0.01);
        All{1,1}{neuron,1}(stimBlock,3) = p;
        All{1,1}{neuron,1}(stimBlock,4) = FSLatency;
        All{1,1}{neuron,1}(stimBlock,5) = jitter;
        All{1,1}{neuron,1}(stimBlock,6) = reliability;
        All{26,1}{neuron,2} = spt_spikeIx;
        try
            [r1,p2] = corrcoef(All{13, 1}(neuron,(45:90)),All{16, 1}{neuron, stimBlock}(1,45:90));
             All{1,1}{neuron,1}(stimBlock,15) = r1(1,2);
        catch
            All{1,1}{neuron,1}(stimBlock,15) = nan;
        end
    end
end

