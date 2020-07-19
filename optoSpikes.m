function [avgSpikes,avgFirstLatency,avgFlashLatency,avgISI,avgAdaptation,avgVarISI] = optoSpikes(filename)

allSpikeNum = [];
allFirstLatency = [];
allFlashLatency = [];
allAvgISI = [];
allAdaptation = [];
allVarISI = [];


spikethreshold = 0; % minimum value for spike
minrepol = -20; % minimum repolarization between spikes
maxwidth = 2; % max spike width

[data,si] = abfload(filename);
for n=1:size(data,3)
    % find data channel based on max variance - only valid if one cell
    [~,dataChannel] = max(var(squeeze(mean(data)),0,2));
    currentSweep = data(:,dataChannel,n);
    flashChannel = size(data,2);
    lightFlashes = findFlashes(data(:,flashChannel,n))*10; % findFlashes decimates
    % find all spikes
    spikebeg = find(currentSweep(2:length(currentSweep)) > spikethreshold & currentSweep(1:length(currentSweep)-1) <= spikethreshold);
   % filter for spikes within light period
    spikebeg = spikebeg(find(spikebeg > lightFlashes(1) & spikebeg < 85000));
%     peaks = zeros(length(spikebeg),2);
%     troughs = zeros(length(spikebeg),2);
    nwidth = maxwidth*1000 / si;
    nspikes = 0;
    spikeLatency = [];
    peaks = [];
    troughs = [];
    for i=1:length(spikebeg) % find spikes, return values & locs of peaks & troughs
        if i < length(spikebeg)
            % find the trough between this spike and the next spike
            [trough, minloc] = min(currentSweep(spikebeg(i):spikebeg(i+1)));
        else
            % find the trough following the final spike in the train
            [trough, minloc] = min(currentSweep(spikebeg(i):spikebeg(i)+2*nwidth));
        end
        % identify the index at which the trough is located
        troughloc = spikebeg(i) + minloc - 1;
      % detemine if the trough is below the minimum repolarization
        % level to detect spikes
        if trough < minrepol
            % if so, determine where the peak of this spike is located
            [peak, peakloc] = max(currentSweep(spikebeg(i):spikebeg(i)+nwidth));
            peakloc = peakloc + spikebeg(i) - 1;
            peaks(i,:) = [peak, peakloc];
            troughs(i,:) = [trough, troughloc];
            flashDiff = lightFlashes - peakloc;
            % find difference between peak and closest preceding lightFlash
            spikeLatency = [spikeLatency abs(max(flashDiff(flashDiff<0)))*si/1000];
        else
            continue;
        end
        nspikes = nspikes + 1;
    end
    
    if isempty(peaks)
        firstLatency = NaN;
        flashLatencies = NaN;
        peaktimes = NaN;
    else
        firstDiff = lightFlashes - peaks(1,2);
        firstLatency = abs(max(firstDiff(firstDiff<0)))*si/1000;
        flashLatencies = mean(spikeLatency);
        peaktimes = peaks(:,2)*si/1000; 
    end
    
    ISIs = diff(peaktimes);

    if isempty(ISIs)
        avgISI = NaN;
    else
        avgISI = mean(ISIs);
    end

    if length(ISIs) > 1
        adaptation = 1/(length(ISIs)-1)*sum((diff(ISIs)./(ISIs(1:end-1) + ISIs(2:end))));
        varISI = var(ISIs);
    else
        adaptation = NaN;
        varISI = NaN;
    end
    
    allSpikeNum = [allSpikeNum nspikes];
    allFirstLatency = [allFirstLatency firstLatency];
    allFlashLatency = [allFlashLatency flashLatencies];
    allAvgISI = [allAvgISI avgISI];
    allAdaptation = [allAdaptation adaptation];
    allVarISI = [allVarISI varISI];
end
    
avgSpikes = [];
avgFirstLatency = [];
avgFlashLatency = [];
avgISI = [];
avgAdaptation = [];
avgVarISI = [];

filters = {[1 0 0 0],[0 1 0 0],[0 0 1 0],[0 0 0 1]};
for i = 1:4
    filter = logical(repmat(filters{i},1,5));
    avgSpikes = [avgSpikes mean(allSpikeNum(filter))];
    avgFirstLatency = [avgFirstLatency mean(allFirstLatency(filter))];
    avgFlashLatency = [avgFlashLatency mean(allFlashLatency(filter))];
    avgISI = [avgISI mean(allAvgISI(filter))];
    avgAdaptation = [avgAdaptation mean(allAdaptation(filter))];
    avgVarISI = [avgVarISI mean(allVarISI(filter))];
end
end
