function [ Sensor_TF, Odata, TFdims, TFtimesMs ] = AGF_Perform_TF( indata, TF, Fdisp, Trials, Etotal, Odata, AnalysisTW)
%%
%performe time-frequency decomposition of the EEG data before running the SVM
%this replaces the time-domain data completely
UseEEGLAB=0;

%Quick check of input
if isfield(TF, 'frequencies')
    if length(TF.frequencies) ~= 2
        disp('Field TF.frequencies must contain a vector of length 2 e.g. [3 25].');return;
    end;
    if ~isfield(TF, 'bandfreq') | isempty(TF.bandfreq)
        disp('TF.bandfreq is missing, substituting with overall frequencies');
        TF.bandfreq = {TF.frequencies};
    end;  
else
    disp('Field TF.frequencies must contain a vector of length 2 e.g. [3 25].');return;
end;
if isfield(TF, 'stepnumber')
    if length(TF.stepnumber) ~= 1
        disp('Field TF.stepnumber must be a scalar of the number of steps, e.g. [20].');return;
    end;
else
    disp('Field TF.stepnumber must be a scalar of the number of steps, e.g. [20].');return;
end;       
if ~isfield(TF, 'smoothing')
    TF.smoothing = 0;
end;     
if ~isfield(TF, 'bands')
    TF.bands = {};
end;     
TWstring='';
if ~isfield(TF, 'time')
    PlotTime = [indata.times(1):indata.times(end)];
    TF.time = PlotTime;
    TWstring = ['Using epoch from ' num2str(PlotTime(1)) ' to ' num2str(PlotTime(end)) ' for TF decomposition - ' num2str(AnalysisTW(1)) ' - ' num2str(AnalysisTW(end)) ' ms will be used for regression.'];
else
    PlotTime = TF.time;
end;   
if ~isfield(TF,'transform') || strcmp(TF.transform, 'none')
    TF.transform = 'none';STrBase = 'Using raw power data';
else
    switch lower(TF.transform)
        case 'normalize'; STrBase = 'Normalizing data within each frequency';
        case 'log'; STrBase = 'Using log-transformed data';
    end;
end;
if ~isfield(TF,'basetype') || strcmp(TF.basetype, 'none')
    TF.basetype = 'none'; STrBase = [STrBase ' and performing no baseline correction.'];
else
    switch lower(TF.basetype)
        case 'percent'; STrBase = [STrBase ' and subtracting single-trial baselines followed by normalization on the average.'];
        case 'subtract'; STrBase = [STrBase ' and subtracting single-trial baselines.'];
    end;
end;
%% 
TFdims  = size(TF.bandfreq,2);

%find time-points with border at edges
if ~isempty(find(indata.times == PlotTime(1)))
    TFt1 = find(indata.times == PlotTime(1));
else
    if Fdisp > 0; warning('Warning: Datalimits are larger than provided data, will use minimum of input data instead...'); end;
    TFt1 = 1; 
end;
if ~isempty(find(indata.times == PlotTime(end)))
    TFt2 = find(indata.times == PlotTime(end));
else
    if Fdisp > 0; warning('Warning: Datalimits are larger than provided data, will use maximum of input data instead...'); end;
    TFt2 = size(indata.times,2);
end;   
TFtimes   = TFt1:TFt2;
Borders   = [find(indata.times == PlotTime(1))-TFt1 TFt2-find(indata.times == PlotTime(end))];%this much border has been added
TFtimesMs = indata.times(TFtimes);
if ~isempty(TF.basetime)
    BaseInd = find(TFtimesMs == TF.basetime(1)):find(TFtimesMs == TF.basetime(2));
    if isempty(BaseInd)
        disp('Error: Baseline appears not to be within the overall time window. Extend the plottime setting (input argument TimeW).');return;
    end;
end;

%%
switch(TF.space) % define if log or line spaced frequencies
    case 'log'
        frex = logspace(log10(TF.frequencies(1)),log10(TF.frequencies(2)),TF.stepnumber);
    case {'linear', 'lin'}
        frex = linspace(TF.frequencies(1),TF.frequencies(2),TF.stepnumber);
end;

%depending on chosen space and step size, not all values will be present: find closest match
for CM = 1 : TFdims
    [~,CloseInd1]       = min(abs(frex-TF.bandfreq{CM}(1)));
    TF.usedfreq{CM}(1)  = frex(CloseInd1);
    [~,CloseInd2]       = min(abs(frex-TF.bandfreq{CM}(2)));
    TF.usedfreq{CM}(2)  = frex(CloseInd2);
    CollapseIndex(CM,:) = [CloseInd1 CloseInd2];
end;
%%
% the width of wavelets scales with the frequency (FWHM of gaussian)
sncy = TF.cyclenumber./(2*pi.*frex);

% make wavelet
t = (indata.xmin-indata.xmax)/2 : 1/indata.srate : abs(indata.xmin-indata.xmax)/2;

for fi=1:length(frex)
    w(fi,:)=exp(2*1i*pi*frex(fi).*t).*exp(-t.^2./(2*sncy(fi)^2)); % sin(2*pi*f*t) IN Euler's formula (e^ik) * gaussian [(-t^2  / SD^s) *2]
end % plot(real(w(1,:)))

Sensor_TF       = nan(length(TFtimes),length(Trials),size(Etotal,1), TF.stepnumber);  %4D time trial channel frequency
Dummy_TF        = nan(length(TFtimes),length(Trials),1,              TF.stepnumber);  %Dummy variable for each electrode to save ram

STrSmooth = '';
if TF.smoothing>0
    Tbef = floor(TF.smoothing/2);
    Tafter = ceil(TF.smoothing/2);
    STrSmooth=['Applying smoothing of ' num2str(Tbef) ' DP before and ' num2str(Tafter) ' DP after each TF time point.'];
    Dummy_TF_bu = Dummy_TF;
end;

                
%%
if Fdisp == 2
    disp(' ')
    disp('********************TF SETTINGS********************')
    disp('*******                                     *******')
    disp(['Frequency range from ' num2str(TF.frequencies(1)) ' to ' num2str(TF.frequencies(2)) ' Hz with ' num2str(TF.stepnumber) ' ' TF.space ' steps.'])
    fprintf([STrSmooth '\n']);
    fprintf([STrBase '\n']);
    fprintf([TWstring '\n\n']);
end;

%%
for chn = 1:length(Etotal)
    chnInd = Etotal(chn);       
    dataX = squeeze(indata.data(chnInd,TFtimes,Trials));
    if UseEEGLAB==1
        [ERSP_Arnaud,~,powbase,Atimes,freqs1,~,~,tfdata] = newtimef( dataX, size(indata.data,2), [TFtimesMs(1) TFtimesMs(end)], indata.srate, TF.cyclenumber, 'freqs', TF.frequencies, 'nfreqs', TF.stepnumber,...
            'freqscale', TF.space, 'baseline', NaN, 'plotitc', 'off', 'plotersp', 'off','scale','abs','trialbase', 'on', 'timesout', 250,'basenorm','off');
        Sensor_TF(:,:,:) = permute(abs(tfdata).^2, [2 3 1]);
        for fi = 1:size(tfdata,1)
            xbase                   = mean(Sensor_TF(1:find(Atimes==-100),:,chn,fi),1);
            Sensor_TF(:,:,chn,fi)  = (Sensor_TF(:,:,chn,fi)-repmat(xbase,size(Sensor_TF,1),1)) ./ mean(xbase);
        end;
    else
        if Fdisp > 0; fprintf('%s\t%s\t','Time-frequency transformation for: ',indata.chanlocs(chnInd).labels); end;
        for fi = 1:TF.stepnumber
            if Fdisp > 0; if rem(fi,5) == 0; fprintf('.'); end; end
            prepdata=fconv_tg(reshape(dataX,1,size(dataX,1)*size(dataX,2)),w(fi,:));    % convolve data with wavelet
            prepdata=prepdata(size(w,2)/2:end-size(w,2)/2);                   % cut of 1/2 the length of the w from beg, and 1/2 from the end
            prepdata=reshape(prepdata,size(dataX,1),size(dataX,2));                     % reshape
            Dummy_TF(:,:,1,fi) = abs(prepdata).^2;                                      % Standard power
            %Use baseline correction for time-frequency data
            switch lower(TF.transform)
                case 'normalize'
                    Dummy_TF(:,:,1,fi)       = normalise(Dummy_TF(:,:,1,fi));
                    STrBase = 'Normalizing data within each frequency';
                case 'log'
                    Dummy_TF(:,:,1,fi)       = log(Dummy_TF(:,:,1,fi));
                    STrBase = 'Using log-transformed data';
                otherwise
                    STrBase = 'Using raw power data';
            end;
            switch lower(TF.basetype)
                case 'percent'
                    xbase                    = mean(Dummy_TF(BaseInd,:,1,fi),1);
                    Dummy_TF(:,:,1,fi)       = (Dummy_TF(:,:,1,fi)-repmat(xbase,size(Dummy_TF,1),1)) ./ mean(xbase); %subtract baseline then divide
                    STrBase = [STrBase ' and subtracting single-trial baselines followed by normalization on the average.'];
                case 'subtract'
                    xbase                    = mean(Dummy_TF(BaseInd,:,1,fi),1);
                    Dummy_TF(:,:,1,fi)       = Dummy_TF(:,:,1,fi)-repmat(xbase,size(Dummy_TF,1),1);
                    STrBase = [STrBase ' and subtracting single-trial baselines.'];
                otherwise
                    STrBase = [STrBase ' and performing no baseline correction.'];
            end; 
            %smooth data if so desired
            if TF.smoothing>0
                Dummy_TF_bu = Dummy_TF; %avoid that smoothing overwrites other datapoints
                for bc = Tbef+1 : length(TFtimes)-Tafter
                    Dummy_TF(bc,:,:,fi) = mean(Dummy_TF_bu(bc-Tbef:bc+Tafter-1,:,:,fi),1);
                end;
            end;
        end;
        Sensor_TF(:,:,chn,:) = Dummy_TF(:,:,1,:);
        if Fdisp > 0; fprintf('\n'); end;
    end;
end

%%
% TFtimesMs=Atimes;
%reduce data over frequency bands (if selected) 
if ~isempty(TF.bands)
    if Fdisp > 0; disp(['Collapsing over ' num2str(TFdims) ' frequency bands...']);end;
    for F =  1 : size(TF.bands,2)
        Sensor_TF_temp(:,:,:,F) = mean(Sensor_TF(:,:,:,CollapseIndex(F,1):CollapseIndex(F,2)),4);
    end;
    Sensor_TF = Sensor_TF_temp; 
    clear Sensor_TF_temp;
else
    TFdims = TF.stepnumber; 
end;

% ExtractTime = find(TFtimesMs==AnalysisTW(1)) : find(TFtimesMs==AnalysisTW(end));
% TFtimesMs = TFtimesMs(ExtractTime);
% Sensor_TF = Sensor_TF(ExtractTime,:,:,:);

Odata.TF            = TF;
Odata.ExtractedTime = TFtimesMs;
Odata.TF.Borders    = Borders;
Odata.TF.freq       = frex;

return;
        
        
        
        
        
        
        
        
        
