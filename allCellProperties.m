% calculating basic cell type properties from a CCIV 
% for standard CCIV - 16 sweeps, 50 pA steps, onset - 200 ms, length - 250
% sweep 5 = 0 current injection

% analyze subthreshold responses @ -50pA (sweep 4)
% analyze individual spikes @ first spike/rheobase
% analyze train responses at rheboase + 1 sweep

% filenames: list of CCIV files as cell array
% eg filenames = {'19410001', '19403001', '19403026', '19403039', '19403052'};
% output file - filename to save output array

function [data,badRows] = allCellProperties(filenames,outputfile)
header = {'Filename', 'Slope', 'AvgFreq', 'MaxFreq', 'Rheobase', 'Tau', 'Rm', 'Vm', 'Halfwidth', 'AHP', 'Peak', 'Upstroke', 'Downstroke', 'rateRatio', 'APthresh', 'Latency', 'firstISI', 'avgISI', 'varISI', 'Adaptation'};
% outputfile = strcat('output/',output);
f = fopen(outputfile,'w');
fprintf(f,'%-15s\t',header{1:end-1});
fprintf(f,'%-15s\n',header{end});

data = cell(length(filenames),length(header));
% badFiles = [];
badRows = [];
for n=1:length(filenames)
    %Establish size, filename
    %current_data = strcat('tempdata/',filenames{n},'.abf')
    data{n,1} = filenames{n};

    % standard inputs for functions
    filename = strcat(filenames{n},'.abf')
    
    try
        [currentData,samplingRate] = abfload(filename);
    catch
        'File not found'
%         badFiles = [badFiles;filename];
        badRows = [badRows; n];
%         data{n,:} = NaN;
        n = n + 1;
%         break
    end
    
    % Checks to make sure file has 16 sweeps for standard CCIV
    if ~ismember(16, size(currentData))
        'Incorrect dimensions, file skipped'
    else

        starttime = 200; % start of current step
        endtime = 450; % end of current step
        zero_sweep = 5; % 0 current injection sweep
        hyper_sweep = 4; % -50 pA current sweep
         spikethreshold = 0; % minimum value for spike
        minrepol = -20; % minimum repolarization between spikes
        maxwidth = 5; % minimum time between spikes
        pulseLength = endtime - starttime;
        currentInjections = -200:50:550;


        % channel with actual data will have larger value than other noise channels 
        [~,dataChannel] = max(abs(mean(currentData)));
        
        % find all spikes 
        [spikeCounts, peaks, troughs] = findSpikes(currentData,samplingRate,dataChannel,spikethreshold,minrepol,maxwidth,starttime,endtime);
        
        try
        % FI curve - get AP train & rheobase measures for future runs
            [FIslope, avgfreq, maxfreq, rheobaseSweep, rheobase_current] = FIcurve(spikeCounts,currentInjections,pulseLength);
            data{n,2} = FIslope;
            data{n,3} = avgfreq;
            data{n,4} = maxfreq;
            data{n,5} = rheobase_current;

             % Hyperpolarized step features
            [tau, Rm, Vrest] = membraneProps(currentData, hyper_sweep, dataChannel, samplingRate, starttime, endtime); 
            data{n,6} = tau;
            data{n,7} = Rm;
            data{n,8} = Vrest;

            % First spike properties
            [halfwidth, ahp, peakheight, peakupstroke, peakdownstroke, rateratio, APthresh] = spikeProperties(currentData, dataChannel, samplingRate, rheobaseSweep, maxwidth, peaks, troughs);     
            data{n,9} = halfwidth;
            data{n,10} = ahp;
            data{n,11} = peakheight;
            data{n,12} = peakupstroke;
            data{n,13} = peakdownstroke;
            data{n,14} = rateratio;
            data{n,15} = APthresh;

            % Spike train properties
            spiking_sweep = rheobaseSweep + 1;
            [latency, firstISI, avgISI, varISI, adaptation] = spikeTrainProperties(peaks, spiking_sweep, starttime, samplingRate);
            data{n,16} = latency;
            data{n,17} = firstISI;
            data{n,18} = avgISI;
            data{n,19} = varISI;
            data{n,20} = adaptation;
        catch
            'Error in spike analysis'
            badRows = [badRows; n];
            n = n + 1;
        end
    end
end



for m = 1:size(data,1) % write data to output file
    fprintf(f,'%-15s\t',data{m,1});
    fprintf(f,'%-15f\t',data{m,2:end-1});
    fprintf(f,'%-15f\n',data{m,end});
end

fclose(f);
end