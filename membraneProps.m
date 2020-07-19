function [tau, Rm, Vrest] = membraneProps(data, sweep, dataChannel,si, starttime, endtime) 
    
currentSweep = data(:,dataChannel,sweep);

% looking at basic membrane properties, smoothing & filtering to get rid of noise
% downsample by a factor of 10
data2 = decimate(currentSweep, 10);

% low-pass filter below 500 Hz
sampfreq = 1 / (si*1e-6);
Wn = 2 * 500 / sampfreq;
data3 = filtfilt(fir1(50,Wn), 1, data2);

%  determine the baseline membrane potential
startindex = starttime / (si*1e-3*10);
baseline = mean(data3(startindex-50:startindex-40));

%  determine the peak hyperpolarization
endindex = endtime / (si*1e-3*10);
peakhyperpol = min(data3(startindex:endindex));

% determine the time constant
decay = min(find(data3(startindex:endindex) < baseline + (peakhyperpol - baseline)*0.632));
tau = decay * si * 10 * 1e-3;

% determine the input resistance
baseline2 = mean(data3(endindex-20:endindex-10));

%Verify that you are using the correct amount of current
    %If using 1st sweep (250pA) divide by 0.25 
    %If using 5th sweep (50pA) divide by 0.05
Rm = (baseline-baseline2) / 0.05; 

Vrest = baseline;