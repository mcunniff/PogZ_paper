function times = findFlashes(sweep)
sweep = round(decimate(sweep,10));
times = [];
for n = 1:length(sweep)-1
    if sweep(n) >=5 && sweep(n+1) <5
        times = [times n];
    end
end

spikeCounts = cell(4,1);
spikeLatency = cell(4,1);


