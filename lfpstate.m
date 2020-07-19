
function [pval_tot, data_info, combdata, coeff2, icasig, incenop2, incencl2, behscore, SCORE, COEFF] = lfpstate(input, behavior)

% what we would like to have is a time series of:
% 1. log power in different frequency bands for each electrode
% 2. amplitude covariation between each frequency band for each
% electrode pair
% 3. phase locking between each frequency band for each electrode
% pair
% 4. cross-frequency coupling between specific frequency bands of
% interest for each electrode pair

% then you would run PCA on the resulting z-scored dataset and
% identify significant dimensions

% the next step is a little unclear.  you could run ICA to identify
% independent components or use k-means clustering to identify
% clusters -- if the latter, you might want to seed using the avg
% activity in the open vs. closed arms, or as mice enter the center
% on runs which avoid vs. explore the open arms, as they exit the
% center, etc.
% -- ?? by seeding does this mean basing the initial centroids on these
% components??

% assuming 'input' is an Nelect x Ntimepts matrix with sampling
% interval dt (in msec)

%% Variable Definition
N = size(input);

dt = 0.0005;

Nelect = N(1);
Ntimepts = N(2);

freqbands = [4 12;
             13 30;
             30 55;
             65 100];

N = size(freqbands);
nbands = N(1);

% CFCslowband = [2 6;
%                6 10];
CFCslowband = [4 12];

CFCfastbands = [13 30;
                30 55;
                65 100];

N = size(CFCfastbands);
ncfc = N(1);

%% Filtering and transforming singals
% first filter and hilbert transform all signals in each frequency band

for i=1:nbands %for each band
    for j=1:Nelect %for each electrode
        % the filter length should be 3 times the period corresponding to
        % the minimum frequency
        L = ceil(3*(1/freqbands(i,1))/dt+1); %finding filter length

        % now filter
        B = fir1(L, [freqbands(i,1) freqbands(i,2)]*dt*2); %filtering

        filtsig{i,j} = filtfilt(B,1,input(j,:));
        x = hilbert(filtsig{i,j}); %hilbert transform
        amp{i,j} = abs(x); %amplitude of the wave
        phase{i,j} = angle(x); %phase of the wave
        
        % amp & phase; have this information for four freq bands over all
        % electrodes

        logamp{i,j} = log(amp{i,j});
    end
end

% how long, in seconds, to use for calculating power, coherence, etc. 
window = 2.5;
win = round(window / dt);

interval = 1;
int = round(interval / dt);

wavfmin = 13;
wavfmax = 100;

n=0;
for f=wavfmin:2:wavfmax
    % compute the Morlet wavelets
    n=n+1;
    sigmat = 1/(2*pi*f/7);
    N = round(2*sigmat / dt);
    t = dt*(-N:N);
    wav{n} = exp(-t.*t / (2*sigmat*sigmat)) .* exp(2*sqrt(-1)*pi*f*t) ...
             / sqrt(sigmat*sqrt(pi));
end

for j=1:Nelect
    fmin = 2;
    fmax = 200;
    L = ceil(3*(1/fmin)/dt+1);        
    B = fir1(L, [fmin fmax]*dt*2);
    y = filtfilt(B,1,input(j,:));
    
    % convolve the boradband filtered signal with wavelets to
    % compute the time varying power at different frequencies
    for n=1:length(wav)
        powertmp = abs(conv(wav{n}, y));
        L = (length(wav{n})-1)/2;
        power{j,n} = powertmp(1+L:length(powertmp)-L);
    end
end

%% Cross-frequency Coherence

N = size(CFCslowband);
    
for k=1:N(1) % cycle through the slow frequency bands for
              % computing CFC
    for i=1:Nelect        
        fmin = CFCslowband(k,1);
        fmax = CFCslowband(k,2);
        L = ceil(3*(1/fmin)/dt+1);        
        B = fir1(L, [fmin fmax]*dt*2);
        
        x = filtfilt(B,1,input(i,:));
        cfcphase{k,i} = angle(hilbert(x));
    end
end

m = 0;

for n=1:win:Ntimepts-win % cycling through all timepoints based on the selected window
    flag = 0;
    for i=1:Nelect
        if ~mean(input(i,n:n+win-1))
            flag = 1;
            break;
        end
    end
    if flag
        break;
    end
    m = m+1;

    % start and end time of this window
    starttime(m) = (n-1)*dt;
    endtime(m) = (n+win-1)*dt;
    
    % compute the mean log power in this window
    for i=1:nbands
        for j=1:Nelect
            % accross the 4 bands, both electrode (4 x 2 matrix of mean log
            % power)
            meanlogamp{i,j}(m) = mean(logamp{i,j}(n:n+win-1));
        end
    end
    
    % for each electrode pair, compute the amp covariation and the wpli
    for i=1:nbands % each band
        l = 0;
        for j=1:Nelect
            for k=j+1:Nelect
                l = l+1;
                
                fmid = mean(freqbands(i,:));
                period = round((1/fmid)/dt);
                rp2 = round(period/2);
                
                x = amp{i,j}(n:n+win-1);
                x = x - mean(x);
                y = amp{i,k}(n:n+win-1);
                y = y - mean(y);
                
                Ctmp = xcorr(x,y,rp2);
                % covariance between waves of the same band between the two
                % electrodes (4 x 1)
                [ampcov{i,l}(m), lag] = max(Ctmp / sqrt(sum(x.*x)*sum(y.*y)));
                
                
                tmp = exp(sqrt(-1)*(phase{i,j}(n:n+win-1)-phase{i,k}(n:n+win-1)));
                X = amp{i,j}(n:n+win-1).*amp{i,k}(n:n+win-1).*tmp;
                imX = imag(X);
                % phase lag index between waves of the same band between
                % the two electrodes (4x1)
                wpli{i,l}(m) = abs(mean(abs(imX) .* sign(imX))) / mean(abs(imX));
            end
        end
    end

    
    % for each electrode pair and freq bands of interest, compute
    % cross-frequency phase-amplitude coupling

    N = size(CFCslowband);
    
    for k=1:N(1)% cycle through the slow frequency bands for
                  % computing CFC
        
        l = 0; % indexes the electrode pair
        f = wavfmin:2:wavfmax;
        for i=1:Nelect
            for j=1:Nelect
                l=l+1;

                % compute the pls for each frequency
                for kk=1:length(wav)
                    plstmp{l}(kk,m) = abs(sum(power{j,kk}(n:n+win) .* exp(-sqrt(-1)*cfcphase{k,i}(n:n+win))) / sum(power{j,kk}(n:n+win)));
                end

                % now combine frequency-specific information into
                % frequency bands
                for kk=1:ncfc % size of fast frequency
                    in = find(f > CFCfastbands(kk,1) & f <= CFCfastbands(kk,2));
                    avgpls{(k-1)*ncfc+kk,l}(m) = mean(plstmp{l}(in,m));
                end
            end
        end
    end
end

%% Creating matrix for PCA
% now we need to z-score (standardize) all these different types of data, and
% then combine them into a single matrix to perform PCA

combdata = [];

n = 0;
for i=1:nbands
    for j=1:Nelect
        n = n+1;
        tmp = meanlogamp{i,j};
        combdata(n,:) = (tmp-mean(tmp))/ sqrt(var(tmp));
        datatype(n) = 0; % provide labels to make it easier to
                         % figure out what different weights in PCs
                         % or ICs represent; datatype = 0
                         % corresponds to power
        datafreq(n) = i;
        datael(n) = j;
    end
end

N = size(ampcov); 
for i=1:N(1)
    for j=1:N(2)
        n=n+1;
        tmp = ampcov{i,j};
        combdata(n,:) = (tmp-mean(tmp))/sqrt(var(tmp));
        datatype(n) = 1; % data type 1 corresponds to amplitude covariance 
        datafreq(n) = i;
        datael(n) = 0;
    end
end

for i=1:N(1)
    for j=1:N(2)
        n=n+1;
        tmp = wpli{i,j};
        combdata(n,:) = (tmp-mean(tmp))/sqrt(var(tmp));
        datatype(n) = 2; % weighted Phase Lag Index
        datafreq(n) = i;
        datael(n) = 0;
    end
end

N = size(avgpls);
for i=1:N(1)
    for j=1:N(2)
        n=n+1;
        tmp = avgpls{i,j};
        combdata(n,:) = (tmp-mean(tmp))/sqrt(var(tmp));
        datatype(n) = 3; % average Phase-Locking Statistic
        datafreq(n) = i+4;
        datael(n) = j+2;
    end
end

%% PCA and finding significant PCs
% now perform PCA on the combined data
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(combdata');

% calculate the Marchenko-Pastur threshold for significant eigenvalues
N = size(combdata);

MP = var(combdata(1:N(1)*N(2))) * (1 + sqrt(N(1)/ N(2)))^2;
sigpcs = sum(LATENT > MP);

%% Correlate PCs with behavioral states
% now determine whether any sig PCs correlate with the two
% behavioral states

% assume 'behavior' is a two row matrix, with row 1 = timepoints,
% and row 2 = zone (e.g. 0 = closed, 1 = center, 2 = open)

% make sure the behavioral vector and data end at the same time
maxtime = min(max(endtime), max(behavior(1,:)));

SCOREin = find(endtime <= maxtime);
SCORE = SCORE(SCOREin,:);

starttime = starttime(SCOREin);

endtime = endtime(SCOREin);

behin = find(behavior(1,:) <= maxtime);
behavior = behavior(:,behin);

n = 0;

% generate a vector of behavior measured in the same time windows
% used before
for i=1:length(SCORE)
    in = find(behavior(1,:) > starttime(i) & behavior(1,:) < ...
              endtime(i));
    if isempty(in) % catches any breaks in the data and assumes same as the previous time frame
        behscore(i) = behscore(i-1);
    else
        % determine whether most timepoints during this window occur in
        % the closed arms vs. center or open arms
        % 0 is closed, 1 is center OR open
        % used to be: behscore(i) = round(mean(behavior(2,in) >= 1));
        behscore(i) = round(mean(behavior(2,in)));
    end
    % behscore --> behavior simplified to be same size as SCORE (ie. both
    % need to be over the same time scale)
end


% compare the correlation between the activity of a PC (or IC) and
% the behavior vector, to the correlation for shuffled data to get
% p-values -- do this instead of just doing a t-test, because
% points in the open (or closed) arms tend to cluster together
% in time

% compare correlation between actual data and random data for p-values to
% avoid issues with data being clusters naturally based on location

R = round(rand(200,1)*length(behscore));

N = size(behscore);

for i=1:sigpcs
    corrval(i) = sum(SCORE(:,i) .* behscore');
    for n=1:200
        shufval(n) = sum(SCORE(:,i)' .* behscore([1+R(n):N(2), 1:R(n)]));
    end
    % p-value = average difference between the actual correlation vs.
    % shuffled/random values, vector with p-values for each sig PC based on
    % where the animal was
    pval(i) = mean(abs(corrval(i)) < abs(shufval));
end

% look for closed->center entries which explore or avoid the open
% arms

mincltime = 5; % minimum durtation (in sec) that a mouse must spend in
               % the closed arms to define this as a closed arm epoch
               % for the purpose of identifying subsequent center
               % entries

minoptime = 0; % minimum duration (in sec) that a mouse must spend
               % in the open arms to define this as a open arm
               % entry (minoptime = 0 just means the mouse has to
               % spend only 1 timepoint in the open arms for this
               % to be the case)

% find all center entries preceded by time spent in the closed arms
incen = 1+find(behavior(2,2:length(behavior)) == 1 & behavior(2,1: ...
                                                  length(behavior)-1) == ...
                                                  0);
% incen --> closed to center, regardless of time spent in closed
x = []; % marks those where it was closed to center, sufficient time in closed
for i=1:length(incen)
    % find the duration of each closed arm epoch
    in2 = max([1,find(behavior(2,1:(incen(i)-1)) > 0)]); % index of when it switches to not closed
    x(i) = behavior(1,incen(i))-behavior(1,in2) > mincltime;
end

% limit to center entries which were preceded by a minimum amount
% of closed arm time
incen = incen(find(x));

% find periods of open arm exploration lasting at least 'minoptime'
% which for open arm is 0, so any open arm exploration
% first find open arm entries 
inop = 1+find(behavior(2,2:length(behavior)) == 2 & behavior(2,1: ...
                                                  length(behavior)- ...
                                                  1) == 1);

x = zeros(1,length(incen));

for i=1:length(inop)
    in2 = min(find(behavior(2,inop(i)+1:length(behavior)) == 1));
    if in2*dt > minoptime
        % find the preceding center entry
        x(max(find(incen < inop(i)))) = 1;
    end
end

incenop = incen(find(x)); % center to open
incencl = incen(find(~x)); % center to closed

% now identify sig PCs which diff during center entries depending on
% whether mice subsequently enter the closed vs. open arms

% find the points in the PCs that correlate not just to position but to the actual act of
% exploring (ie. center to closed vs. center to open)
incenop2 = [];
incencl2 = [];
for i=1:length(incenop)
    % find the index within the PC (or IC) which corresponds to the
    % index of the behavior vector  
    if behavior(1, incenop(i)) <= starttime(length(starttime))
        incenop2(i) = min(find(starttime > behavior(1,incenop(i))));
    else
        incenop2(i) = length(starttime);
    end
end

for i=1:length(incencl)
    % find the index within the PC (or IC) which corresponds to the
    % index of the behavior vector
    if behavior(1,incencl(i)) <= starttime(length(starttime))
        incencl2(i) = min(find(starttime > behavior(1,incencl(i))));
    else
        incencl2(i) = length(starttime);
    end
end

% compute p-values that significant PCs differ for center entries
% that subsequently explore vs. avoid the open arms

for i=1:sigpcs
    [h,pval2(i)] = ttest2(SCORE(incenop2,i), SCORE(incencl2,i));
end

%% Calculate ICs, run same stats on ICs
% project onto the significant PCs; find ICs; run statistics on ICs

[icasig, A, W] = fastica(SCORE(:,1:sigpcs)'); 

% just running ICA on the scores of the significant PCs
% means that ICA is pulling out the independent components of the most
% important components of the data, defined by running PCA

sigics = size(icasig,1); % only different when one of the ics doesn't converge

% this is finding the p-value for the correlation between IC and the behavior
for i=1:sigics % used to be sigpcs --> but sometimes IC doesn't converge
    corrval2(i) = sum(icasig(i,:) .* behscore);
    for n=1:200
        shufval(n) = sum(icasig(i,:) .* behscore([1+R(n):N(2), 1:R(n)]));
    end
    pval3(i) = mean(abs(corrval2(i)) < abs(shufval));
end


% compute p-values that ICs differ for center entries
% that subsequently explore vs. avoid the open arms

for i=1:size(icasig,1)
    [h,pval4(i)] = ttest2(icasig(i,incenop2), icasig(i,incencl2));
end

pval_tot = cell(1,4);
pval_tot{1} = pval;
pval_tot{2} = pval2;
pval_tot{3} = pval3;
pval_tot{4} = pval4;

% compute the coefficients of the ICs in terms of the original data
% have found the weight matrix for the ICs -- basically the shift from what
% is really happening to the signal observed 
coeff2 = W*COEFF(:,1:sigpcs)';

% could cluster the ICs across mice to identify corresponding ICs in
% different mice; pool statistics for ICs in the same clusters in
% different mice

%% Feature Detection

% makes sig_feat_ind where each row is a different PC, column is the index
% of the feature in the original data, ranked by coefficients 

% makes sig_feat_type where each row is a different PC, column is the datatype of
% the feature, ranked by coefficients
% where: 0 --> power
        % 1 --> amplitude covariance
        % 2 --> phase lag index (quantifying phase synchronization)
        % 3 --> phase-locking statistic (quantifying cross-frequency
        % phase-amplitude coupling)
        
% makes sig_feat_freq where each row is a different PC, column is the
% frequencies studied in that feature, ranked by coefficients
% where: 1 --> [4 12]
        % 2 --> [13 30]
        % 3 --> [30 55]
        % 4 --> [65 100]
        % 5 --> [13 30], [2 6]; beta and theta
        % 6 --> [30 55], [2 6]; low gamma and theta
        % 7 --> [65 100], [2 6]; high gamma and theta
        % 8 --> [13 30], [6 10]; beta and alpha
        % 9 --> [30 55], [6 10]; low gamma and alpha
        % 10 --> [65 100], [6 10]; high gamma and alpha
        
% makes sig_feat_el where each row is a different PC, column is the
% electrode for that feature, ranked by coefficients
% where 0 --> between HPC and mPFC
        % 1 --> HPC
        % 2 --> mPFC
        % 3 --> HPC to HPC
        % 4 --> HPC to mPFC
        % 5 --> mPFC to HPC
        % 6 --> mPFC to mPFC
            

coeff2_temp = abs(coeff2);
for i = 1:size(icasig,1)
    for j = 1:length(coeff2)
        curr_ind = find(coeff2_temp(i,:) == max(coeff2_temp(i,:)));
        sig_feat_ind(i,j) = curr_ind;
        sig_feat_type(i,j) = datatype(curr_ind);
        sig_feat_freq(i,j) = datafreq(curr_ind);
        sig_feat_el(i,j) = datael(curr_ind);
        coeff2_temp(i,curr_ind) = -1;
    end
end

data_info = cat(3, datatype, datafreq, datael);
sorted_data_info = cat(3, sig_feat_type, sig_feat_freq, sig_feat_el, sig_feat_ind);


end

