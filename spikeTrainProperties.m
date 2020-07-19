function [latency, firstISI, avgISI, varISI, adaptation] = spikeTrainProperties(allPeaks, sweep, starttime, si)

peaks = allPeaks{sweep,1};

% latency from beginning of current pulse to first spike
latency = (peaks(1,2)*si/1000 - starttime)/1000;
peaktimes = peaks(:,2)*si/1000; 
% find interspike interval - distance between all spikes
ISIs = diff(peaktimes);

if isempty(ISIs)
    firstISI = NaN;
    avgISI = NaN;
    varISI = NaN;
else
    % interval between first two spikes
    firstISI = ISIs(1);
    % average interspike interval
    avgISI = mean(ISIs);
    % variance in interspike interval
    varISI = std(ISIs)/mean(ISIs);
end

if length(ISIs) > 1
    % adaptation - whether spikes get further apart as train goes on
    adaptation = 1/(length(ISIs)-1)*sum((diff(ISIs)./(ISIs(1:end-1) + ISIs(2:end))));
else
    adaptation = NaN;
end

%  
%  
end
 
 
 