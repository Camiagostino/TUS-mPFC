function [ FreqOut, OutBase ] = AGF_Hilbert( indata, baseline, frequencies, tx_in, start, ende, srate, basemethod)
%%
%Function that performs Hilbert Transformation on EEG data.
%Thanks Jim Cavanagh.
%
%Input:
%'indata'       EEG data or icaacct structure (channels, datapoints, epochs)
%'baseline'     Baseline used (milliseconds [-1000 -500])
%'frequencies'  Frequencies to be analyzed
%'tx'           Timepoints in the data
%'srate'        Sampling rate
%'basemethod'   subtract = just subtract baseline value of each trial
%               divide = subtract and then divide each trial
%
%
%Output:
%FreqOut        4 Dimensional awway with:
%               1 = Frequency band, 2 = channel, 3 = datapoints, 4 = epochs
%               

dims=size(indata);
%srate=500; % sampled at 500 Hz...
tx=tx_in(1):2:tx_in(2);  % ...so data are -1500 to 1500 ms peri-feedback
t1=find(tx==start); % time range start window
t2=find(tx==ende); % time range end window
%clear FreqOut Power Power_base
% Split into different bands, also get raw EEG
for freq = 1 : length(frequencies)
    if strcmp(lower(frequencies{freq}),'delta'),     lo=1; hi=4;    % 1 = Delta
    elseif strcmp(lower(frequencies{freq}),'theta'), lo=4; hi=8;    % 2 = Theta
    elseif strcmp(lower(frequencies{freq}),'alpha'), lo=8; hi=12;   % 3 = Alpha
    elseif strcmp(lower(frequencies{freq}),'beta'), lo=12; hi=20;   % 4 = Beta
    elseif strcmp(lower(frequencies{freq}),'gamma'), lo=20; hi=50;  % 5 = Gamma
    elseif strcmp(lower(frequencies{freq}),'raw'), lo=[]; hi=20;    % 6 = Raw EEG
    end

    % Best to filter in two stages: first low pass (high band)
    temp1=AGF_eegfilt(indata,srate,[],hi);
    if freq<6
       % If not raw EEG, high pass filter (low band) as well
       temp2=AGF_eegfilt(temp1,srate,lo,[]);
       % Get power envelope:
       temp3=abs(hilbert(temp2'))'.^2;  % HILBERT OPERATES COLUMN-WISE.  REQUIRES TRANSPOSITION AND BACK!
    else
       % for raw EEG, no power envelope
       temp3=temp1;
    end
    Power=reshape(temp3,dims(1),dims(2),dims(3));
    clear Power_base
    if size(baseline,1)==1
        b1=find(tx==baseline(1)); % baseline start window
        b2=find(tx==baseline(2)); % baseline end window
    end;
    for trial = 1 : dims(3)
        if size(baseline,1)>1
            if mod(baseline(trial,1),2)
                b1=find(tx==baseline(trial,1)-1);
            else
                b1=find(tx==baseline(trial,1));
            end;
            if mod(baseline(trial,2),2)
                b2=find(tx==baseline(trial,2)-1);
            else
                b2=find(tx==baseline(trial,2));
            end;
        end;
        
        %Get single trial baseline
        bsl_ST(trial,:)=mean(Power(:,b1:b2,trial),2);
        bsl=repmat(mean(Power(:,b1:b2,trial),2),1,length(t1:t2));
        %Scale by baseline
        Power_base(:,:,trial)=Power(:,t1:t2,trial)-bsl;
    end;
    if strcmp(basemethod,'divide')
        %Divide by average activity in baseline timewindow
        Power_base=Power_base./repmat(abs(mean(bsl_ST))',[1 length(t1:t2) size(Power,3)]);
    end;
    % Now data will be 4D: with the first dim as freq above
    % (the first 5 are band-specific power envelopes, the 6th is filtered EEG)
    % The other 3D are as the EEG.data structure (chans, time, epochs)
    FreqOut(freq,:,:,:)=Power_base;
    clear temp* Pow*;
end;
OutBase = mean(bsl_ST);
return
% close all
% figure
% hold on
% plot([-1000:2:1000],squeeze(mean(FreqOut(2,1,:,:),4)),'g')
% %plot([-1000:2:1000],squeeze(mean(FreqOut(6,19,:,:),4)) ./ 10,'k') % Raw EEG scaled by factor of 10
% set(gca,'xlim',[-1000 1000]);
