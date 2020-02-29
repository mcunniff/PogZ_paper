% Returns power and coherence measures around a single timepoint (eg for
% use in analyzing transitions) 
% Can be expanded to multiple electrodes but currently output optimized for
% 2 electrodes/1 pair, such as PFC/HPC recordings

% Modfied by Margaret Cunniff from ptscript3

% THESE ARE THE PARAMATERS THAT NEED TO BE SET BY THE USER
% filelist, freqband, channels, and uniquetypes

% WT
% filelist = 'dec16_pogZ_WT_EPM_runTypes.txt';
% dataID = 'dec16_wt_runTypes';

% HET
% filelist = 'dec16_pogZ_Het_EPM_runTypes.txt';
% dataID = 'dec16_het_runTypes';

% VIP 
filelist = 'archRealTimeRunTypes.txt';
dataID = 'feb20_VIP_';

% freqband = [4 12; 12 30; 30 58; 62 90];
freqband = [1 4; 4 12];

lofreq = [4 12]; % the low frequency band (in Hz)

hifreq = [62 90]; % the high frequency band (in Hz)

% channels = [2 3]; % which electrodes are to be analyzed

types = [0,1,2,3]; % the types associated with timepoints
               
fid = fopen(filelist);

eegfile = fscanf(fid, '%s', [1]);

% N = size(channels);

% nelpairs = N(1);

N = size(lofreq);
ncfcbands = N(1);

n = 0;

ntypes = length(types);


while eegfile,
    disp('Current File')
    disp(eegfile)
    disp(' ')
    timefile = fscanf(fid, '%s', [1]);
    
    if ~timefile,
        break;
    end
    
    n = n+1;
    nspect(n,:) = zeros(1,ntypes);
    
    [data, header] = readedf5(eegfile);

%     channels = channelIDs(n,:);
    if header.channelname(1,1:3) == 'BIO'
        channels = [2 3];
    else
        channels = [1 2];
    end
    
    if size(data{1,channels(1)}) ~= size(data{1,channels(2)})
        disp('Warning: Data channels are not same dimensions')
    end
%     delt = header.duration / header.nsamples(channels(1)));
    channel1 = data{1, channels(1)};
    channel2 = data{1, channels(2)};

    % check if data are the same size. if not, downsample larger one.
    if length(channel1) == length(channel2)
        delt = header.duration / header.nsamples(channels(1));

    elseif length(channel1) > length(channel2)
        factor = length(channel1)/length(channel2);
        data{1,channels(1)} = downsample(channel1,factor);
        delt = header.duration / header.nsamples(channels(2));

    elseif length(channel2) > length(channel1)
        factor = length(channel2)/length(channel1);
        data{1,channels(2)} = downsample(channel2,factor);
        delt = header.duration / header.nsamples(channels(1));
    end

    if size(data{1, channels(1)}) ~= size(data{1, channels(2)})
        msg = 'Data still not equal, downsampling failed';
        error(msg)
    end

    % cycle through all the timepoints
    fid2 = fopen(timefile);
    timepts = fscanf(fid2, '%f', [3,inf]);
    logical = [timepts(1,:)>7.5; timepts(1,:)>7.5; timepts(1,:)>7.5];
    timepts = reshape(timepts(logical),3,[]);
    N = size(timepts);
    


    % how many points on either side of the timepoint to analyze?
    npts = 7.5 / delt; % examine activity for +/- 7.5 sec
    
    % window size to use for computing spectrogram and coherogram
    win = 1 / delt; % use a 1-second window

    % window size to use for computing pls, amplitude covariation,
    % and cross-frequency coupling
%     winpls = 2.5 / delt; % use a 2.5-second window
     % FOR DELTA
     winpls = 2.5 / delt;
    

    for i=1:N(2),
%         disp(timepts(:,i))
%         
%         % first convert the timepoint (in seconds) to indices
        in = round(timepts(1,i)/delt) - 1;
        type = timepts(3,i) ;
        typein = find(types == type);
        indivRefPoints{n,i} = in;
        if in-npts < 1 || in+npts > length(data{channels(1)}),
            disp('Time outside bounds of data')
            break;
        end
        
      
        % now compute the spectrogram in a window of time around
        % this point
        [S,F,T] = spectrogram(data{channels(1)}(in-npts+1:in+npts),win, ...
                        [],[],1/delt);
        if nspect(n,typein),
            S1mean{n,typein} = S1mean{n,typein} + abs(S);
        else
            S1mean{n,typein} = abs(S);
        end        
        
        S = spectrogram(data{channels(2)}(in-npts+1:in+npts),win, ...
                        [],[],1/delt);
        if nspect(n,typein),
            S2mean{n,typein} = S2mean{n,typein} + abs(S);
        else
            S2mean{n,typein} = abs(S);
        end

        % next compute the coherence around this point
        [C,F2] = mscohere(data{channels(1)}(in-npts+1:in+npts),data{channels(2)}(in-npts+1:in+npts),win,[],[],1/delt);

%         [WC, ~,F3] = wcohere(data{channels(1)}(in-npts+1:in+npts),data{channels(2)}(in-npts+1:in+npts),seconds(1/delt));
        
        if nspect(n,typein),
            Cmean{n,typein} = Cmean{n,typein} + abs(C);
        else
            Cmean{n,typein} = abs(C);
        end

        % now compute the amplitude covariation and phase locking
        % around this point

        N2 = size(freqband);
        for j=1:N2(1),
            l = 0;
            for k=in-npts+1:(winpls/2):in+npts-winpls+1,                
                l = l+1;
                x = data{channels(1)}(round(k):round(k)+winpls-1);
                y = data{channels(2)}(round(k):round(k)+winpls-1);
                [maxval,lagtmp,plstmp,plsphase,wplitmp] = laglead6(x,y,delt,freqband(j,1), freqband(j,2));
                % HERE
                if nspect(n,typein),
                    plsmean{n,typein}(j,l) = plsmean{n,typein}(j,l) + ...
                        plstmp;
                    ampcor{n,typein}(j,l) = ampcor{n,typein}(j,l) + ...
                        maxval;
                    lag{n,typein}(j,l) = lag{n,typein}(j,l) + ...
                       lagtmp;                   
                    wpli{n,typein}(j,l) = wpli{n,typein}(j,l) + ...
                        wplitmp;                    
                else
                    plsmean{n,typein}(j,l) = plstmp;
                    ampcor{n,typein}(j,l) = maxval;
                    lag{n,typein}(j,l) = lagtmp;
                    wpli{n,typein}(j,l) = wplitmp;
                end
                indivPLS{n,typein}(j,l,i) = plstmp;
                indivAmp{n,typein}(j,l,i) = maxval;
                indivLag{n,typein}(j,l,i) = lagtmp;                   
                indivWpli{n,typein}(j,l,i) = wplitmp;  
            end
        end
        % finally compute the cross-frequency coupling around this point

        l = 0;
        for j=1:ncfcbands,
            for k=in-npts+1:(winpls/2):in+npts-winpls+1,                
                l = l+1;
                x = data{channels(1)}(round(k):round(k)+winpls-1); % HPC
                y = data{channels(2)}(round(k):round(k)+winpls-1); % PFC
                
                [C,F2] = mscohere(x,y,win,[],[],1/delt);
% 
% %         [WC, ~,F3] = wcohere(data{channels(1)}(in-npts+1:in+npts),data{channels(2)}(in-npts+1:in+npts),seconds(1/delt));
%         
                if nspect(n,typein),
                    Cmean2{n,typein}(l,:) = Cmean2{n,typein}(l,:) + abs(C)';
                else
                    Cmean2{n,typein}(l,:) = abs(C)';
                end
            
% 
%                 [phasestmp, plstmp, ampprofile] = pac3(x,y,delt,lofreq(j,:),hifreq(j,:));
% 
%                 if nspect(n,typein),
%                     cfcmean1{n,typein}(j,l) = cfcmean1{n,typein}(j,l) ...
%                         + plstmp;
%                     meanampprofile1{n,typein}(l,:) = ...
%                         meanampprofile1{n,typein}(l,:) + ampprofile';
%                 else
%                     cfcmean1{n,typein}(j,l) = plstmp;
%                     meanampprofile1{n,typein}(l,:) = ampprofile';
%                 end
%                 indivCfc1{n,typein}(j,l,i) = abs(mean(plstmp));
%                 indivAmpProfile1{n,typein}{j,l,i} = ampprofile;
%                 indivPhases1{n,typein}{j,l,i} = phasestmp;
                [phasestmp, plstmp, ampprofile] = pac3(y,x,delt,lofreq(j,:),hifreq(j,:));

                if nspect(n,typein),
                    cfcmean2{n,typein}(j,l) = cfcmean2{n,typein}(j,l) + plstmp;
                    meanampprofile2{n,typein}(l,:) = ...
                        meanampprofile2{n,typein}(l,:) + ampprofile';
                else
                    cfcmean2{n,typein}(j,l) = plstmp;
                    meanampprofile2{n,typein}(l,:) = ampprofile';
                end
                indivPhases2{n,typein}{j,l,i} = phasestmp;
                indivCfc2{n,typein}(j,l,i) = abs(mean(plstmp));
                indivAmpProfile2{n,typein}{j,l,i} = ampprofile;
                
                indivTimePts{n,typein}(j,l,i) = k;
            end
        end

        if nspect(n,typein),
            nspect(n,typein) = nspect(n,typein) + 1;
        else,
            nspect(n,typein) = 1;
        end

    end
    
% 
%     for i=1:ntypes,
%         S1mean{n,i} = S1mean{n,i} / nspect(n,i);
%         S2mean{n,i} = S2mean{n,i} / nspect(n,i);
%         Cmean2{n,i} = Cmean2{n,i} / nspect(n,i);
%     end
%     
% % %         cfcmean1{n,i} = cfcmean1{n,i} / nspect(n,i);
% %         cfcmean2{n,i} = cfcmean2{n,i} / nspect(n,i);
% % %         meanampprofile1{n,i} = meanampprofile1{n,i} / nspect(n,i);
% %         meanampprofile2{n,i} = meanampprofile2{n,i} / nspect(n,i);
% %         plsmean{n,i} = plsmean{n,i} / nspect(n,i);
% %         ampcor{n,i} = ampcor{n,i} / nspect(n,i);
% %         lag{n,i} = lag{n,i} / nspect(n,i);
% %         wpli{n,i} = wpli{n,i} / nspect(n,i);
%     end
%     
    eegfile = fscanf(fid, '%s', [1]);
     clear delt header i in j k l lagtmp m npts win winpls
end


clear ans C channels data delt eegfile fid filelist freqband header
clear hifreq i in j k l lagtmp lofreq m n N N2 ncfcbands nelpairs
clear npts nspect ntypes  S timefile timepts type typein
clear types win winpls x y plsphase  wplitemp maxval
clear  fid2 logical plstemp ampprofile phasetmp

t = linspace(-7.5,7.5,11);
save(dataID);