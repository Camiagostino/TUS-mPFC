function [ All_tf, Power_Trial, frex ] = AGF_Tfreq( EEG, elect, ICA, num_freqs, sfreq, efreq, space, numcycles, t1, t2, t3, t4 )
%UNTITLED Theo TF Function
%   Detailed explanation goes here
if exist('fconv_JFC','file') ~=2,
    error('function: fconv_JFC missing');
end

switch(space) % define if log or line spaced frequencies
    case 'log'
        frex = logspace(log10(sfreq),log10(efreq),num_freqs);
    case 'lin'
        frex = linspace(sfreq,efreq,num_freqs);
end

 % the width of wavelets scales with the frequency (FWHM of gaussian)
sncy = numcycles./(2*pi.*frex);

% make wavelet
t = size(EEG.times,2)/2/1000*-1:1/EEG.srate:size(EEG.times,2)/2/1000;

for fi=1:length(frex)
    w(fi,:)=exp(2*1i*pi*frex(fi).*t).*exp(-t.^2./(2*sncy(fi)^2)); % sin(2*pi*f*t) IN Euler's formula (e^ik) * gaussian [(-t^2  / SD^s) *2]
end % note: to view a single wavelet: plot(real(w(1,:)))

% go trough frequencies
T1 = find(EEG.times==t1);
T2 = find(EEG.times==t2);
if ~isempty(t3) & ~isempty(t4)
    T3b = find(EEG.times==t3);
    T4b = find(EEG.times==t4);
    T3=T3b-T1+1;
    T4=T4b-T1+1;
end;
%elect = find(strcmpi({EEG.chanlocs.labels},elect));

All_tf=nan(num_freqs,length(T1:T2),3);
Power_Trial=nan(num_freqs,length(T1:T2),size(EEG.data,3));
if ICA==1
    dataX = squeeze(EEG.icaact(elect,:,:));
else
    dataX = squeeze(EEG.data(elect,:,:));
end;
fprintf(1,['Progress: Percent of frequencies analyzed:   ']);
dims = size(dataX);
for fi=1:num_freqs

    do=floor( fi / num_freqs *100);
    if do < 10
        fprintf(2,'\b%d', do);
    else
        fprintf(2,'\b\b%d', do);
    end;   

    stuffnjunk=fconv_JFC(reshape(dataX,1,dims(1)*dims(2)),w(fi,:));                 % convolve data with wavelet
    stuffnjunk=stuffnjunk((floor((size(w,2)-1)/2):end-1-ceil((size(w,2)-1)/2)));    % cut of 1/2 the length of the w from beg, and 1/2 from the end
    stuffnjunk=reshape(stuffnjunk,dims(1),dims(2));
    Power_Trial(fi,:,:) = abs(stuffnjunk(T1:T2,:)).^2;
    All_tf(fi,:,2) = mean(Power_Trial(fi,:,:),3)';
    %All_tf(fi,:,2) = mean(abs(stuffnjunk(T1:T2,:)).^2,2);                           % Standard power
    All_tf(fi,:,3) = abs(mean(exp(1i*(angle(stuffnjunk(T1:T2,:)))),2));             % Lachaux PLV

    if ~isempty(t3) & ~isempty(t4)
        basevalST   = mean(abs(stuffnjunk(T3b:T4b,:)).^2,1);                            % Get baseline activity in the single trials
        baseval     = mean(All_tf(fi,T3:T4,2),2);

        %Power_Trial(fi,:,:) = 10*log10(bsxfun(@rdivide, squeeze(Power_Trial1), basevalST)); % Single Trial Baseline correction standardized at every trial + dB
        Power_Trial(fi,:,:) = 10*log10(Power_Trial(fi,:,:) ./ baseval);                      % Single Trial Baseline correction standardized at average + dB
        All_tf(fi,:,2)      = 10*log10( All_tf(fi,:,2) ./ baseval);                          % Baseline correction + dB
        Power_Trial(fi,:,:) = Power_Trial(fi,:,:) - mean(mean(Power_Trial(fi,T3b:T4b,:),2),3);
    else
        Power_Trial(fi,:,:) = 10*log10((abs(stuffnjunk(T1:T2,:)).^2));
    end;
%         
%         clc
%         STbasevalDP1 = mean(basevalST)
%         basevalDP1   = baseval
%         STvalueDP1   = mean(abs(stuffnjunk(1,:)).^2)
%         normalDP1    = mean(abs(stuffnjunk(1,:)).^2,2)
%         
%         disp('output')
%         STvalueDP1   = mean(Power_Trial(fi,1,:),3)
%         normalDP1    = All_tf(fi,1,2)
    clear stuffnjunk stuffn2junk;

end
All_tf(1,:,1) = mean(EEG.data(elect,T1:T2,:),3);
fprintf(2,'\n')
end