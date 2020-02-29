function [phases, pls, ampprofile] = pac3(x,y,dt,flo,fhi)

% estimates phase-amplitude coupling in two signals
% x is the signal that will be filtered in the high frequency band
% then the amplitude envelope of x in the high frequency band will
% be filtered in the low frequency band to extract the phase
% then the phase of the hi freq amplitude envelope of x will be
% compared with y in the lo freq band

% phase is a list of all the phase differences
% pls is a list of phase locking values calculated using a fixed
% size time window

% flo is a 2 element array with the upper/lower bound of the low
% frequency band

% fhi is a 2 element array with the upper/lower bound of the high
% frequency band 

% the filter length should be 3 times the period corresponding to
% the minimum frequency
L = ceil(3*(1/fhi(1))/dt+1);
L2 = ceil(3*(1/flo(1))/dt+1);

% now generate the filters
B = fir1(L, [fhi(1) fhi(2)]*dt*2);
B2 = fir1(L2, [flo(1) flo(2)]*dt*2);

% filter the first signal at the higher frequency band

x2 = filtfilt(B,1,x);

% now compute the hilbert transform and extract the amplitude
% envelope of the fast oscillation

x3 = hilbert(x2);

ampx = abs(x3);

% now filter the fast oscillation amplitude signal in the lower frequency band

ampxfilt = filtfilt(B2, 1, ampx);
y2 = filtfilt(B2, 1, y);

% now compute the phases in the lower frequency band

phasex = angle(hilbert(ampxfilt));
phasey = angle(hilbert(y2));

% now compute the correlations as a function of time and lag/lead
N = length(phasex);

% also compute the amplitude of the fast oscillation as a function
% of the phase of the slow oscillation

nbins = 6;

ampprofile = zeros(nbins,1);
ampcounts = zeros(nbins,1);

for i=1:N,
    phi = mod(phasey(i), 2*pi);
    phasebin = ceil(phi / (2*pi/nbins));
    ampprofile(phasebin) = ampprofile(phasebin) + ampx(i);
    ampcounts(phasebin) = ampcounts(phasebin) + 1;
end

ampprofile = ampprofile ./ ampcounts;

% compute the number of timesteps corresponding to the central
% frequency - then throw away the first 3 and last 3 periods worth
% of data since there may be edge effects of filtering in these periods
period = round((1/mean(flo))/dt);

in = (1+3*period):(N-period*3);

% compute the phase differences
phases = phasex(in)-phasey(in);

% the window size for computing the strength phase-locking should
% be 10 times the lowest frequency

win = round(10*(1/mean(flo))/dt);
stepsize = win;

N = length(phases);

pls = [];

% now compute many estimates of phase locking using the fixed
% window size

npls = 0;
for n=1:stepsize:N-win,
    npls = npls+1;
    pls(npls) = mean(exp(2*pi*sqrt(-1)*phases(n:n+win)));
end
