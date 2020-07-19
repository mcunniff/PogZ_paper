function [spikeCounts, peaks, troughs] = findSpikes(data,si,dataChannel,spikethreshold,minrepol,maxwidth,starttime,endtime)
% FILENAME is the abffile to analyze
% THRESHOLD is the sweep that represents the threshold, i.e. the
% first nonzero response
% NSWEEPS is the number of sweeps
% SPIKETHRESHOLD is the voltage level (in mV) for detecting spikes -- not
% the actual spike threshold which will be determined.  Probably
% should be something like -15 mV.
% MINREPOL is a minimum voltage that must be attained between
% spikes to "reset" the spike detector (probably should be ~ -30
% mV)
% MINREFRACTORY is the minimum refractory period (in msec) between detected
% spikes (usually ~2 msec)
% STARTTIME is the time in msec during the sweep to start the
% analysis, e.g. the time at which the current pulse begins
% ENDTIME is the time when the current pulse ends

% step through all sweeps of a recording and find all spikes
numSweeps = size(data,3);
spikeCounts = zeros(numSweeps,1);
peaks = cell(numSweeps,1);
troughs = cell(numSweeps,1);

for k=1:numSweeps
    currentSweep = data(:,dataChannel,k);
    N = length(currentSweep);
    
    % initialize vars for sweep stats
    nspikes = 0;
    currentPeaks = [];
    currentTroughs = [];
    
    % find times at which the membrane potential crosses the voltage threshold for detecting spikes
    spikebeg = find(currentSweep(2:N) > spikethreshold & currentSweep(1:N-1) <= spikethreshold);

    % determine which of these potential spikes occur during the current pulse
    endtime = endtime + 25;
    spikebeg = spikebeg(find(spikebeg*si/1000 > starttime & spikebeg*si/1000 < endtime));

    for i=1:length(spikebeg) % find spikes, return values & locs of peaks & troughs

        if i < length(spikebeg)
            % find the trough between this spike and the next spike
            [trough, minloc] = min(currentSweep(spikebeg(i):spikebeg(i+1)));
        else
            % find the trough following the final spike in the train
            [trough, minloc] = min(currentSweep(spikebeg(i):spikebeg(i)+2*maxwidth*1000/si));
        end

        % identify trough location relative to beginning of sweep
        troughloc = spikebeg(i) + minloc - 1;

        % detemine if the trough is below the minimum repolarization level 
        if trough < minrepol
            [peak, peakloc] = max(currentSweep(spikebeg(i):spikebeg(i)+maxwidth*1000/si));
            % identify peak location relative to beginning of sweep
            peakloc = peakloc + spikebeg(i) - 1;
            currentPeaks = [currentPeaks; peak, peakloc];
            currentTroughs = [currentTroughs; trough troughloc];
        else
            continue;
        end
        nspikes = nspikes + 1;
    end
    spikeCounts(k,1) = nspikes;
    peaks{k,1} = currentPeaks;
    troughs{k,1} = currentTroughs;
end

