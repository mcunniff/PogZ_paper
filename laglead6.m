function [maxval,lag,pls,plsphase,wpli] = laglead6(x,y,dt,fmin,fmax)

% first calculate the central frequency
fmid = 0.5*(fmin+fmax);

% the filter length should be 3 times the period corresponding to
% the minimum frequency
L = ceil(3*(1/fmin)/dt+1);

% now filter
B = fir1(L, [fmin fmax]*dt*2);

x2 = filtfilt(B,1,x);
y2 = filtfilt(B,1,y);

% now compute the hilbert transforms of both filtered signals and
% extract the amplitude envelopes

x3 = hilbert(x2);
y3 = hilbert(y2);

ampx = abs(x3);
ampy = abs(y3);

phasex = angle(x3);
phasey = angle(y3);

ampx2 = ampx - mean(ampx);
ampy2 = ampy - mean(ampy);

period = round((1/fmid)/dt);
rp2 = round(period/2);

m = 0;
for j=-rp2:rp2,
    m = m+1;
    if j < 0,
        acx(m) = sum(ampx2(1:length(ampx2)+j).* ampx2(1:length(ampx2)+j));
        acy(m) = sum(ampy2(-j+1:length(ampy2)).* ampy2(-j+1:length(ampy2)));
    else
        acx(m) = sum(ampx2(j+1:length(ampx2)) .* ampx2(j+1:length(ampx2)));
        acy(m) = sum(ampy2(1:length(ampy2)-j).* ampy2(1:length(ampy2)-j));
    end
end

Ctmp = xcorr(ampx2,ampy2,rp2);

Ctmp = Ctmp ./ sqrt(acx'.*acy');

[maxval,lag] = max(Ctmp);
lag = lag*dt;
tmp = exp(2*pi*sqrt(-1)*(phasex-phasey)); % differences between phase at each time point
plstmp = mean(tmp); % mean of phase difference - normal PLS
plsphase = angle(plstmp); % frequency of PLS
pls = abs(plstmp); % amplitude of PLS

X = ampx.*ampy.*tmp; % multiplying all phase differences by amplitude at that point
imX = imag(X); % taking the imaginary component
wpli = abs(mean(abs(imX) .* sign(imX))) / mean(abs(imX)); % 
