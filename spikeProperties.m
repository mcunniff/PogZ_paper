function [halfwidth, ahp, peakheight, peakupstroke, peakdownstroke, rateratio, APthresh] = spikeProperties(data, dataChannel, si, sweep, maxwidth, allPeaks, allTroughs)

% HALFWIDTH is the time (in msec) for the membrane potential to rise from the
% point halfway between the peak and trough to the peak, and back
% to the halfway point
% AHP is the AHP membrane potential
% PEAKHEIGHT is the height of the spike
% MAXRISE is the maximum rising slope (in mV/timestep)
% MAXFALL is the maximum falling slope (in mV/timestep)
% PEAKTIMES are the times at which each spike occurs (more
% precisely when the peak of each spike occurs)
% THRESH is the threshold, determined as the point at which the
% third derivative of the membrane potential is maximal
% ISIs is a vector with all the interspike intervals
% VREST is the resting membrane potential, just the average
% membrane potential before the current pulse

nwidth = maxwidth*1000/si;
currentSweep = data(:,dataChannel,sweep);

% select peaks & troughs for sweep of interest
peaks = allPeaks{sweep,1};
troughs = allTroughs{sweep,1};

% basic properties of first spike
peakheight = peaks(1,1);
ahp = troughs(1,1);

% halfwidth of first spike - height = distance from peak to trough
peakloc = peaks(1,2);
troughloc = troughs(1,2);

% going backwards from peak, looking for when it goes below half height
for j1=peakloc:-1:peakloc-nwidth
    if currentSweep(j1) < 0.5*(peakheight + ahp)
        break;
    end
end

% going forward from peak, looking for when it goes below half height
for j2=peakloc:peakloc+nwidth 
    if currentSweep(j2) < 0.5*(peakheight + ahp)
        break;
    end
end

halfwidth = ((j2-j1)*si/1000);

% peak upstroke & downstroke - look for maximum amount of change between two points, normalize by dt
% gradient - calculates difference between two adjacent points

spikeROI = currentSweep(peakloc-nwidth:troughloc);
[~,B] = max(gradient(gradient(spikeROI))); % max of second derivative - spike threshold
APthresh = spikeROI(B);
threshloc = find(spikeROI(B) == spikeROI);

riseROI = currentSweep(threshloc+(peakloc-nwidth):peakloc);
fallROI = currentSweep(peakloc:troughloc);

peakupstroke = max(diff(riseROI))/si; % deriv in mV/ms
peakdownstroke = min(diff(fallROI))/si;
rateratio = abs(peakupstroke/peakdownstroke);

end
