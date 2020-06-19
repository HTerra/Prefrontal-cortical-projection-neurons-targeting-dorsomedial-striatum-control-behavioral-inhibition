function spikes = get_spikes_from_kwik(fileName,cluster)
%% get_spikes_from_kwik: return spike timestamps array of kwik file
% Return the spike timestamps in the kwikfile specified by fileName input
% argument and the cluster (ie cell/unit) (use 'all' instead of cluster nr
% to return all spiketimes. Returns spike timestamps in
% seconds in an array.
RecParam;
time_samples = hdf5read(fileName, '/channel_groups/0/spikes/time_samples');
belongsToCluster = hdf5read(fileName, '/channel_groups/0/spikes/clusters/main');
if ischar(cluster) && strcmp(cluster,'all')
   spikeSamples = time_samples; 
else
    spikeSamples = time_samples(belongsToCluster == cluster);
end
spikes = (double(spikeSamples)-1)/samplingRate; %Convert to secs, do -1 because sample 1  is 0 seconds so 30001 is 1 second
spikes = spikes'; % Return a horizontal vector
end
