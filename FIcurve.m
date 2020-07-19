function [slope, avgfreq, maxfreq, rheobaseSweep, rheobase] = FIcurve(spikeCounts, currentInjections, pulseLength)
% spikeCounts is array of spike counts for each sweep
% zero_sweep is # of sweep with no current injection

% find first sweep with spikes
rheobaseSweep = find(spikeCounts,1);
% edge cases where spontaneous spiking causes rheobase to be artifically
% low
if rheobaseSweep < 6
    rheobaseSweep = 6;
end

% find current injection at rheobase - sweeps are in 50pA increments
rheobase = currentInjections(rheobaseSweep);
numSweeps = length(spikeCounts);
above_thresh = spikeCounts(rheobaseSweep:numSweeps);

%find average firing frequency & max firing frequency
avgfreq = (sum(above_thresh)/length(above_thresh))/(pulseLength/1000);
maxfreq = max(spikeCounts)/(pulseLength/1000);

if (numSweeps - rheobaseSweep) <= 4
    slopeSweep = numSweeps;
    slopeSpikes = above_thresh;
else
    slopeSweep = rheobaseSweep + 3;
    slopeSpikes = above_thresh(1:4);
end

% fit line for spike count vs current injection to get FI slope
X = [(currentInjections(rheobaseSweep:slopeSweep))' ones(size(slopeSpikes))];
slopes = regress(slopeSpikes, X);
slope = slopes(1);




