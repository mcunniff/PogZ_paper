function [avgPeaks,avgFirstRatio,avgLastRatio,avgArea,avgDelay,avgDecay,avgShape,delayVar] = optoIPSC(filename,dataChannel)
[data,si] = abfload(filename);
allPeaks = [];
allLocs = [];
firstRatios = [];
lastRatios = [];
areas = [];
peakDelays = [];
decayTimes = []; 
normShape = []; 

for n=1:size(data,3)
    % find data channel based on max variance - only valid if one cell
    if dataChannel == 0
        [~,dataChannel] = max(var(squeeze(mean(data)),0,2));
    end
    currentSweep = data(:,dataChannel,n);
    currentSweep = decimate(currentSweep,10);
    currentSweep = medfilt1(currentSweep,10);
    
    flashChannel = size(data,2);
    lightFlashes = findFlashes(data(:,flashChannel,n));
    numFlashes = length(lightFlashes);
    flashInterval = lightFlashes(2) - lightFlashes(1);

    baseline = mean(currentSweep(1:4000));
    endBase = mean(currentSweep(8001:12000));
    rest = (baseline + endBase)/2;
    
    normalized = (currentSweep - rest);
    auc = abs(trapz(normalized(6001:8000)));
    
    noise = 2.5*std(normalized)+mean(normalized);
    
    [peakVal, peakLoc] = max(normalized);
    
    halfFlashes = numFlashes/2;
    if rem(halfFlashes,2) == 0
    
        allFlashes = zeros(halfFlashes,1);
        allLocs = zeros(halfFlashes,1);
        for k = halfFlashes:numFlashes
            if k < numFlashes
                [allFlashes(k), allLocs(k)] = max(normalized(lightFlashes(k-1+halfFlashes):lightFlashes(k+halfFlashes)));
            else
                [allFlashes(k), allLocs(k)] = max(normalized(lightFlashes(k-1+halfFlashes):lightFlashes(k+halfFlashes)+flashInterval));
            end
            allLocs(k) = allLocs(k) + lightFlashes(k);
        end
        normFlashes = abs(allFlashes./allFlashes(1));
        allNorm = (sum(normFlashes))/numFlashes;    

        if abs(allFlashes(1)) > abs(noise)
            peakDelay = (allLocs(1)-lightFlashes(1))*si/100;

            decayVal = (1/exp(1))*allFlashes(1);
            decayLoc = find(normalized(allLocs(1):length(normalized)-1) > decayVal & normalized(allLocs(1)+1:length(normalized)) <= decayVal);
            if decayLoc(1) < (lightFlashes(2) - lightFlashes(1))
                decayTime = decayLoc(1)*si/100;
            else
                decayTime = NaN;
            end
        %     
            firstRatio = allFlashes(2)/allFlashes(1);
    %         lastRatio = allFlashes(numFlashes)/allFlashes(1);

            normPeak = abs(peakVal);
            allPeaks = [allPeaks normPeak];
            firstRatios = [firstRatios firstRatio];
    %         lastRatios = [lastRatios lastRatio];
            areas = [areas auc];
            peakDelays = [peakDelays peakDelay];
            decayTimes = [decayTimes decayTime];
            normShape = [normShape allNorm];
        else
            allPeaks = [allPeaks 0];
            firstRatios = [firstRatios NaN];
    %         lastRatios = [lastRatios NaN];
            areas = [areas 0];
            peakDelays = [peakDelays NaN];
            decayTimes = [decayTimes NaN];
            normShape = [normShape NaN];
        end
    else
        allPeaks = [allPeaks 0];
        firstRatios = [firstRatios NaN];
%         lastRatios = [lastRatios NaN];
        areas = [areas 0];
        peakDelays = [peakDelays NaN];
        decayTimes = [decayTimes NaN];
        normShape = [normShape NaN];
    end
end
avgPeaks = [];
avgFirstRatio = [];
avgLastRatio = [];
avgArea = [];
avgDelay = [];
avgDecay = [];
avgShape = [];
delayVar = [];

filters = {[1 0 0 0],[0 1 0 0],[0 0 1 0],[0 0 0 1]};
for i = 1:4
    filter = logical(repmat(filters{i},1,5));
    avgPeaks = [avgPeaks mean(allPeaks(filter))];
    avgFirstRatio = [avgFirstRatio mean(firstRatios(filter))];
%     avgLastRatio = [avgLastRatio mean(lastRatios(filter))];
    avgArea = [avgArea mean(areas(filter))];
    avgDelay = [avgDelay mean(peakDelays(filter))];
    avgDecay = [avgDecay mean(decayTimes(filter))];
    avgShape = [avgShape mean(normShape(filter))];
    delayVar = [delayVar var(peakDelays(filter))];
end