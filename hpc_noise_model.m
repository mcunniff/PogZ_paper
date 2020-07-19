% timestep in msec
dt = 0.1;

hpcfreq = 0.5; % theta oscillations frequency in Hz
hpcrate = 100; % Peak poisson rate of hippocampal spikes
noiserate = 50;

% membrane time constants in msec
taupyr = 10;
tauint = 10;

% simulation duration in msec
T = 1000;

% excitatory synaptic weights, for hippocampal spikes or noise spikes onto
% each cell type
hpcpyr = 0.6;
hpcint = 0.33;
noisepyr = 0.6;
noiseint = 0.33;

% inhibitory synaptic weights
intpyr = 0.1;
vipint = 0.1;

% synaptic time constants in msec
pyrexctau = 8;
pyrinhtau = 20;

intexctau = 8;
intinhtau = 20;

vipexctau = 8;

hpcWeights = 0.1:0.05:1;
weights = hpcWeights;
noiseWeights = 1;
numIterations = 1000;

finalHpcCorrMean = zeros(length(hpcWeights),length(noiseWeights));
finalNoiseCorrMean = zeros(length(hpcWeights),length(noiseWeights));
finalPyrSpikesMean = zeros(length(hpcWeights),length(noiseWeights));
finalFsSpikesMean = zeros(length(hpcWeights),length(noiseWeights));
finalRatioMean = zeros(length(hpcWeights),length(noiseWeights));

finalHpcCorrStd = zeros(length(hpcWeights),length(noiseWeights));
finalNoiseCorrStd = zeros(length(hpcWeights),length(noiseWeights));
finalPyrSpikesStd = zeros(length(hpcWeights),length(noiseWeights));
finalFsSpikesStd = zeros(length(hpcWeights),length(noiseWeights));
finalRatioStd = zeros(length(hpcWeights),length(noiseWeights));

finalHpcCorrVals = cell(length(hpcWeights),length(noiseWeights));
finalNoiseCorrVals = cell(length(hpcWeights),length(noiseWeights));
finalPyrSpikeCounts = cell(length(hpcWeights),length(noiseWeights));
finalFsSpikeCounts = cell(length(hpcWeights),length(noiseWeights));


for j = 1:length(hpcWeights)
    hpcint = hpcWeights(j);
    disp(['HPC Weight: ',num2str(hpcint)])
    for k = 1
        if k == 1
            noiseint = hpcint;
        else
            
            noiseint = 0.33;
        end
        disp(['Noise Weight: ', num2str(noiseint), ', Trial: 1'])
        allHpcCorr = [];
        allNoiseCorr = [];
        allPyrSpikes = [];
        allFsSpikes = [];
        for i = 1:numIterations
            if rem(i,100) == 0
                disp(['Weight: ', num2str(hpcint), ', Trial: ', num2str(i)])
            end
            % synaptic variables (these are updated at each timestep and indicate how
            % active each synapse is
            pyrexcsyn = 0;
            pyrinhsyn = 0;
            intexcsyn = 0;
            intinhsyn = 0;
            vipexcsyn = 0;

            % variables representing membrane potentials for each cell type
            pyrvm = 0;
            intvm = 0;
            vipvm = 0;

            % keep track how many spikes of each type have occurred
            npyrspikes = 0;
            nintspikes = 0;
            nvipspikes = 0;
            nhpcspikes = 0;
            nnoisespikes = 0;

            % keep track of the times of spikes of each type
            pyrtimes = [];
            inttimes = [];
            viptimes = [];
            hpctimes = [];
            noisetimes = [];

            % thresholds for each cell type, in mV
            pyrthresh = 10;
            intthresh = 10;
            vipthresh = 10;

            % calculate how many timesteps there will be
            N = ceil(T/dt) + 1;

            % random numbers used to calculate when hippocampal or noise spikes occur
            R = rand(N,2);

            t = 0:dt:T;

            % assume that the hippocampal spike rate oscillates based on theta
            % frequency
            theta = 0.5*sin(t*2*pi*hpcfreq/1000)+0.5;

            n=0; % counter for how many timesteps have elapsed
            
            for t=0:dt:T,
                n=n+1;
                % voltages undergo passive decay
                pyrvm = pyrvm - pyrvm*dt/taupyr;
                if pyrvm < -10
                    pyrvm = -10;
                end
                intvm = intvm - intvm*dt/tauint;
                vipvm = vipvm - vipvm*dt/tauvip;

                % synaptic variables undergo passive decay
                pyrexcsyn = pyrexcsyn - pyrexcsyn*dt/pyrexctau;
                pyrinhsyn = pyrinhsyn - pyrinhsyn*dt/pyrinhtau;
                intexcsyn = intexcsyn - intexcsyn*dt/intexctau;
                intinhsyn = intinhsyn - intinhsyn*dt/intinhtau;
                vipexcsyn = vipexcsyn - vipexcsyn*dt/vipexctau;

                % voltage changes driven by synaptic variables
                pyrvm = pyrvm + hpcpyr*pyrexcsyn - intpyr*pyrinhsyn;
                intvm = intvm + hpcint*intexcsyn - vipint*intinhsyn;
                vipvm = vipvm + hpcvip*vipexcsyn;

                % determine if there is a hippocampal spike and update excitatory synapses
                if R(n,1) < hpcrate*theta(n)*dt/1000,
                    nhpcspikes = nhpcspikes + 1;
                    hpctimes(nhpcspikes) = t;
                    pyrexcsyn = pyrexcsyn + hpcpyr;
                    intexcsyn = intexcsyn + hpcint;
                    vipexcsyn = vipexcsyn + hpcvip;
                end

                % determine if there is a noise spike and update excitatory synapses    
                if R(n,2) < noiserate*dt/1000,
                    nnoisespikes = nnoisespikes + 1;
                    noisetimes(nnoisespikes) = t;
                    pyrexcsyn = pyrexcsyn + noisepyr;
                    intexcsyn = intexcsyn + noiseint;
                    vipexcsyn = vipexcsyn + noisevip;
                end

                % determine if the neurons spike and update synapses
                if pyrvm > pyrthresh,
                    if npyrspikes == 0
                        npyrspikes = npyrspikes + 1;
                        pyrtimes(npyrspikes) = t;
                        pyrvm = 0;
%                         pyrTrace = [pyrTrace 40];
%                         pyrRaster = [pyrRaster 1];
                    elseif pyrtimes(npyrspikes) + 15 < t
                        npyrspikes = npyrspikes + 1;
                        pyrtimes(npyrspikes) = t;
                        pyrvm = 0;
                    end
                end

                if intvm > intthresh,
                    if nintspikes == 0
                        nintspikes = nintspikes + 1;
                        inttimes(nintspikes) = t;
                        pyrinhsyn = pyrinhsyn + 1;
                        intvm = 0;

                    elseif inttimes(nintspikes) + 5 <= t
                        nintspikes = nintspikes + 1;
                        inttimes(nintspikes) = t;
                        pyrinhsyn = pyrinhsyn + 1;
                        intvm = 0;

                    end
                end

            end


            x = 0:2.5:1000;

            % create a filter to smooth the spike trains
            filterwidth = 10; % width of the filter in msec
            x2 = -2*filterwidth:2.5:filterwidth;
            y2 = exp(-x2.*x2/(2*filterwidth*filterwidth));
            y2 = y2 / sum(y2); % normalize the filter so its elements sum to 1

            % bin the spikes
            yhpc = histc(hpctimes, 0:2.5:1000);
            ypyr = histc(pyrtimes, 0:2.5:1000);
            ynoise = histc(noisetimes, 0:2.5:1000);

            % smooth the binned spike trains using the filter
            yhpc2 = filtfilt(y2,1,yhpc);
            ypyr2 = filtfilt(y2,1,ypyr);
            ynoise2 = filtfilt(y2,1,ynoise);
            % plot(x, xhpc2, x, xpyr2)

            % compute the correlations between the smoothed spike trains
            if ~isempty(yhpc2) && ~isempty(ypyr2) && ~isempty(ynoise2)
                hpcCorr = corr(yhpc2', ypyr2');
                noiseCorr = corr(ynoise2', ypyr2');
                allHpcCorr = [allHpcCorr hpcCorr];
                allNoiseCorr = [allNoiseCorr noiseCorr];
            else
                allHpcCorr = [allHpcCorr NaN];
                allNoiseCorr = [allNoiseCorr NaN];
            end

            allPyrSpikes = [allPyrSpikes npyrspikes];
            allFsSpikes = [allFsSpikes nintspikes];
    
        end
        allRatios = allHpcCorr./allNoiseCorr;
        finalRatioMean(j,k) = mean(allRatios);
        finalRatioStd(j,k) = std(allRatios);

        finalHpcCorrVals{j,k} = allHpcCorr;
        finalHpcCorrMean(j,k) = nanmean(allHpcCorr);
        finalHpcCorrStd(j,k) = std(allHpcCorr);

        finalNoiseCorrVals{j,k} = allNoiseCorr;
        finalNoiseCorrMean(j,k) = nanmean(allNoiseCorr);
        finalNoiseCorrStd(j,k) = std(allNoiseCorr);

        finalPyrSpikeCounts{j,k} = allPyrSpikes;
        finalPyrSpikesMean(j,k) = mean(allPyrSpikes);
        finalPyrSpikesStd(j,k) = std(allPyrSpikes);

        finalFsSpikeCounts{j,k} = allFsSpikes;
        finalFsSpikesMean(j,k) = mean(allFsSpikes);
        finalFsSpikesStd(j,k) = std(allFsSpikes);

    end
end
finalRatio = finalHpcCorrMean./finalNoiseCorrMean;
clearvars -except weights final*