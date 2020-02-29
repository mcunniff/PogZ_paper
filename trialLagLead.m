function [covar,covarlag,pls,plsphase,wpli] = trialLagLead(x,y,dt,fmin,fmax)
% Calculating covariance, PLS, and WPLI for a given trial
% maxval - max amplitude covariation
% lag - lag at which max amplitude covariation occurs
% pls - magnitude of phase locking statistic
% plsphase - phase of PLS
% wpli - weighted phase locking index

% Modified by Margaret Cunniff from laglead6 (Sept 2017)
% x - data from electrode 1
% y - data from electrode 2
% dt - sampliing period of data
% fmin - lower bound of frequency window
% fmax - upper bound of frequency window


% first calculate the central frequency
fmid = 0.5*(fmin+fmax);

% filter order - here 3 times the length of the minimum frequency
% order of window - number of samples used as an input to filter
L = ceil(3*(1/fmin)/dt+1);

% fir1 creates a hamming filter window
% creates bandpass filter between fmin and fmax, calculated over L samples
B = fir1(L, [fmin fmax]*dt*2);

% filtering the 2 signals using the hamming window
% B: filter, 1: specfies B as FIR (filter type)
x2 = filtfilt(B,1,x);
y2 = filtfilt(B,1,y);

% computing hilbert transform of filtered data - used to calculate envelope
% and extract instantaneous phase and amplitude from the signal
x3 = hilbert(x2);
y3 = hilbert(y2);

% instantaneous amplitude of signal
ampx = abs(x3);
ampy = abs(y3);

% instantaneous phase of signal
phasex = angle(x3);
phasey = angle(y3);

% normalizing instantaneous amplitudes
normAmpX = ampx - mean(ampx);
normAmpY = ampy - mean(ampy);

    
% Comparing signals shifted to each other -  rp2 +/- is max offset (so half
% a cycle in either direction based on middle of freq band)
period = round((1/fmid)/dt);
rp2 = round(period/2);

% take the non-normalized cross correlation between x & y for delays of +/-
% rp2
rawCorr = xcorr(normAmpX,normAmpY,rp2);

% eg for a given time i and delay the corr would be X[i+delay]*Y[i]
% that would be summed over all possible values of i such that vectors are the same length:
%(eg X(delay:len(X)) & Y(1:len(Y)-delay)
% to normalize, use autocorrelation for the data included in each
% correlation calculation (so accounting for delay changing total length)
lenX = length(normAmpX);
lenY = length(normAmpY);
allLags = [-rp2:rp2];
for shift=1:length(allLags)
    lag = allLags(shift); 
    if lag >= 0        
        autoX(shift) = sum(normAmpX(lag+1:lenX).^2);
        autoY(shift) = sum(normAmpY(1:lenY-lag).^2);
    else
        autoX(shift) = sum(normAmpX(1:lenX+lag).^2);
        autoY(shift) = sum(normAmpY(-lag+1:lenY).^2);
    end
end

% normalizing cross correlation by the auto correlation of the data
% corresponding to that delay


normCorr = rawCorr ./sqrt(autoX'.*autoY');

[covar,covarlag] = max(normCorr);
% maxval = maximum cross correlation = amplitude covariation number
% lag = position at which maxval occurs = phase offset with max covariation
covarlag = covarlag*dt;

% convert the phase difference at each time point into a complex number
% with an amplitude of 1 and angle equal to the phase diff
phaseDiffs = exp(sqrt(-1)*(phasex-phasey)); 
% take mean of all these vectors, left with single complex number for mean PLS 
plsVector = mean(phaseDiffs);
% take the angle of this vector for the phase of PLS and magnitude of the
% vector for the PLS value
plsphase = angle(plsVector);
pls = abs(plsVector);

% multiple instantaneous amplitudes and difference in instantaneous phase
crossSpectra = ampx.*ampy.*phaseDiffs; 
 % taking the imaginary component - in complex plane, y value
imagCrossSpectra = imag(crossSpectra);
%Larger phase difference -> larger imaginary component
% Small phase difference -> sensitive to noise
% WPLI weights those with larger imaginary compnents   
% absolute value imaginary number = magnitude 
% multiple the magnitude of the imaginary component by the sign,
wpli = abs(mean(abs(imagCrossSpectra) .* sign(imagCrossSpectra))) / mean(abs(imagCrossSpectra)); % 
