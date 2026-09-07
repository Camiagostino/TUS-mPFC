function [ALL] = quantify_ERPs_singletrial(inname, inlocation, MinMax, zeit, electrode, VP, varargin)
%Function that calculates ERPs and returns each value.
%
%Required:
%inname: Name of inputfile
%inlocation: Location of inputfile
%MinMax:     Calculate Mean, Minimum or Maximum in first timewindow
%zeit:       Timewindow(x:y as datapoints)
%electrode:  Number of Electrode to calculate with. ATTENTION: If more than
%            one Electrode is specified, function will first cluster electrodes by
%            calculating mean over all specified electrodes, and then perform normal
%            calculation on the result. Input has to be a Vector like [13 22 30].
%VP:         Subjects as array of strings {'01';'02'...}
%
%Optional: 
%If >1 intervals are defined, timewindows and way of calculation to process both have to be defined.
%'MinMax2'  = [min,max,mean,minmean,maxmean]
%'zeit2'    = [x:y]
%'action'   = ['a'/'s'/'m'/'n'] a=add, s=subtrakt, m=mean, n=nothing performed in order of 'order'.
%'relative' = [0/1/2] Activates relative Mode. Timespan 1 is used as absolute
%             timepoint from where onwards other times may be used. If set to '2', the
%             next timepoint will be relative to the previously found one, if '1' all are relative to the first one. 
%             For Example: To calculate the FRN you can set time1 to 200:250 and mode to 'max'. If relative
%             is activated and time2 is 50 and mode2 is 'min', the function returns the
%             minimal value of the next 50 timepoints behind the found first peak. If
%             relative is set to 2 and a third timewindow is specified, a third peak can
%             be searched in the third timewindow, which begins after the end of the second.
%             Function will always calculate mean of absolute and last relative peak and use 'aktion' to
%             calculate, all peaks are returned with their relative times.
%             Relative timepoints can be positive or negative!
%'TimeMode' = [X] default is 1: Times are returned in relation
%             to first point of dataset. If X is a number, this is assumed as the
%             timepoint from wich onward time is calculated in ms assuming 500Hz
%             recording and used as 'X0' for plotting.
%'Hz'       = [X] can be set to change frequency to calculate time as ms
%             (default = 500) from TimeMode onwards.
%'plot'     = [0,1] if set to 1, a plot indicating all found peaks is produced but not autamaticly saved
%'EqualYLim'= [0,1] sets the YLimits for the plot in all subjects to the same value (1), or uses automatic 
%             values (0).
%'InclEpoch'= Numer of epochs to include. Must be a cell array with the same size as VP.
%'OrdEpoch' = ['rnd' / 'dis'] Order of epochs if less Epochs than all are specified. 'rnd' = random sample of epchs, 'dis' = equal distribution of
%             all epochs.
%'RetSTD'   = Return std of whole epoch (1), of certain timeperiod [X:Y] or nothing (0, default)
%'RetTrnr'  = Default(0) does not return, (1) returns EEG.epoch.trnr for each measurement (useful if you want to combine behavioural data).
%'Surround' = Only necessary if measurement is set to minmean or maxmean. Describes datapoints around found peak to calculate mean on.

%AG Fischer, 2010

%SET DEFAULT VALUES
nargs = nargin-6;
inputNum=1;         %Standard is 1
Frelative=0;    	%Default setting for relative mode
FTimeMode=1;        %Calculte time in ms from datapoint 1
FHz=500;            %Dafault 500Hz
Fplot=0;            %Do not plot by default
order= [1 2];       %Default order is 1:3
action=['n' 'n'];   %Default action is to do nothing (n)
res1=NaN;           %Result of Calculation step 1
res2=NaN;           %Result of Calculation step 2
FEqualYLim=0;       %Automatic YLimits for plotting.
ESET=0;             %Include all epochs.
FRetSTD=0;          %Do not return standard deviation.
FRetTrnr=0;         %Do not return Trialnumbers.
Fsurround=0;        %No surrounding timewindow should be used.

if nargs > 1 
  if ~(round(nargs/2) == nargs/2)
    error('Odd number of input arguments??')
  end
end;

for i = 1:2:length(varargin)
    Param = varargin{i};
    if ~isstr(Param)
      error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
     case 'aktion'
         Faktion=varargin{i+1};
     case 'minmax2'
         FMinMax2=varargin{i+1};inputNum=2;      
     case 'minmax3'
         FMinMax3=varargin{i+1};inputNum=3;
     case 'zeit3'
         Fzeit3=varargin{i+1}; 
     case 'zeit2'
         Fzeit2=varargin{i+1};
     case 'relative'
         Frelative=varargin{i+1};
     case 'plot'
         Fplot=varargin{i+1}; 
     case 'action'
         action=varargin{i+1};
     case 'inclepoch'
         Finclepoch=varargin{i+1};
         ESET=1;
     case 'ordepoch'
         Fordepoch=varargin{i+1}; 
     case 'retstd'
         FRetSTD=varargin{i+1}; 
     case 'rettrnr'
         FRetTrnr=varargin{i+1}; 
     case 'surround'
         Fsurround=varargin{i+1}; 
     case 'timemode'
         FTimeMode=varargin{i+1};
         if FTimeMode<1
             display('*Starting timepoint must be positive!')
             return;
         end;
     case 'hz'
         FHz=varargin{i+1};
         if isstr(FHz)
            display('*Hz must be a number!')
            return;
         end;
     case 'equalylim'
         FEqualYLim=varargin{i+1};
    end;
end; %Get Variables from varargin array

%LOOP THROUGH VP
for vpl=1:length(VP)
    clear EEG Abweichung time1 amp1 time1 L time1TP time1ms TimeWindow2 amp2 time2 tim2TP time2ms res1 all_amps all_times EpochSample
    res1=NaN;
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset( 'filename',['pb' VP{vpl} inname], 'filepath', inlocation);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = eeg_checkset( EEG );
    
    %Reduce number of Epochs if only a sample should be taken
    if ESET==1        
        if Finclepoch{vpl}>length(EEG.epoch)
            display('*WARNING: Input number to reduce epochs is longer than number of epochs.')
        elseif strcmp(Fordepoch,'rnd') 
            display('*Using random sample of epochs.')
            Epochsample=randperm(Finclepoch{vpl});
        elseif strcmp(Fordepoch,'dis')
            display('*Using independent sample of epochs.')
            fak=length(EEG.epoch)/Finclepoch{vpl};
            Epochsample=1:Finclepoch{vpl};
            Epochsample=floor(Epochsample(:)*fak);
        end;
        EEG.data=EEG.data(:,:,Epochsample);
        EEG.epoch=EEG.epoch(Epochsample); %Reduce Epochs too (in order to return Trnr correctly)
    end;
    
    %Cluster Electrodes if more than one electrode is specified
    if length(electrode)>1
        EEG.data(electrode(1),:,:)=mean(EEG.data(electrode,:,:));
        display(['*Electrodes are clustered over: ' num2str(electrode)])
        Felectrode=electrode(1);clear electrode;electrode=Felectrode;
    end;

    %If set, return std for Epoch or period in epoch
    if FRetSTD~=0
        if size(FRetSTD,2)>1
            display(['*Returning STD of each epoch from datapoint: ' num2str(FRetSTD(1)) ' until ' num2str(FRetSTD(end))]);
            Abweichung=squeeze(std(EEG.data(electrode,FRetSTD,:)))';
        else
            display('*Returning STD of each epoch.');
            Abweichung=squeeze(std(EEG.data(electrode,:,:)))';
        end;
    else
        Abweichung=NaN;
    end;
    
    %ALWAYS PERFORM CALCULATION OF FIRST PEAK IN zeit1
    if strcmp(MinMax,'min')
        [amp1,time1]=min(EEG.data(electrode,zeit,:));display('*Taking Minimum in first timewindow')
    elseif strcmp(MinMax,'max')
        [amp1,time1]=max(EEG.data(electrode,zeit,:));display('*Taking Maximum in first timewindow')
    elseif strcmp(MinMax,'mean')
        amp1=mean(EEG.data(electrode,zeit,:));display('*Taking mean in first timewindow') 
        time1=NaN; 
    elseif strcmp(MinMax,'minmean')
        if Fsurround<1;display('*Error, surrounding times must be positive');return;end;
        [temp1,time1]=min(EEG.data(electrode,zeit,:));
        ttemp1=squeeze(((time1+min(zeit)-1)));
        amp1=mean(EEG.data(electrode,[ttemp1(1)-Fsurround:ttemp1(1)+Fsurround],:));
        display('*Searching minimum and taking surrounding times mean.');clear ttemp1 temp1
    elseif strcmp(MinMax,'maxmean')
        if Fsurround<1;display('*Error, surrounding times must be positive');return;end;
        [temp1,time1]=min(EEG.data(electrode,zeit,:));
        ttemp1=squeeze(((time1+min(zeit)-1)));
        amp1=mean(EEG.data(electrode,[ttemp1(1)-Fsurround:ttemp1(1)+Fsurround],:));
        display('*Searching maximum and taking surrounding times mean.');clear ttemp1 temp1 
    end;

    %Transform to absolute values in dataset and milliseconds
    time1TP=squeeze(((time1+min(zeit))-1))';
    time1ms=squeeze(((time1+min(zeit))-FTimeMode-1)*(1000/FHz))';
    amp1=squeeze(amp1)';
    L=length(time1TP); %Used just to shorten coding

    %Define second timewindow and then calculate and transform second peak
    if inputNum>=2 
        if Frelative>0 %Timewindow2 relative to first peak
            if length(Fzeit2)>1;display('*Error, relative times cannot be intervals');return;end;
            if strcmp(MinMax,'mean');display('*Error, absolute timewindow 1 (zeit) must be set to "min" or "max" if relative.');return;end;
            display('*Using relative timespans! First set Timespan is used as absolute...')
            for c=1:L
                if Fzeit2>0
                    TimeWindow2(c,:)=time1TP(c)+1:time1TP(c)+Fzeit2;
                else
                    TimeWindow2(c,:)=time1TP(c)+Fzeit2:time1TP(c)-1;
                end;
            end;
        else %Timewindow2 has absolute value
            if length(Fzeit2)<2;display('*Error, absolute times must be intervals >2');return;end;
            for c=1:L
                TimeWindow2(c,:)=Fzeit2;
            end;
        end;

        %Calculate second peak
        for c=1:L
            if strcmp(FMinMax2,'min')
                [amp2(c),time2(c)]=min(EEG.data(electrode,TimeWindow2(c,:),c));
            elseif strcmp(FMinMax2,'max')
                [amp2(c),time2(c)]=max(EEG.data(electrode,TimeWindow2(c,:),c)); 
            elseif strcmp(FMinMax2,'mean')
                [amp2(c)]=mean(EEG.data(electrode,TimeWindow2(c,:),c));
                time2(c)=NaN; 
            end;    
        %Transform to absolute values in dataset and milliseconds
        time2TP(c)=squeeze(((time2(c)+min(TimeWindow2(c,:)))-1));
        time2ms(c)=squeeze(((time2(c)+min(TimeWindow2(c,:)))-FTimeMode-1)*(1000/FHz));
        end;
        amp2=squeeze(amp2);
    end;

    %Perform calculation if defined
    %Check if enough actions for the datasets are specified
    if length(order)<length(action)-1 | length(action)<inputNum-1
        display('Error: Not enough actions specified for order and input datasets!')
        return;
    end;

    if action(1)=='a' & inputNum>1
        res1=eval(['amp' num2str(order(1))])+eval(['amp' num2str(order(2))]);
        display('Adding first two peaks...')
    elseif action(1)=='s'
        res1=eval(['amp' num2str(order(1))])-eval(['amp' num2str(order(2))]);
        display('Subtracting first two peaks...')
    elseif action(1)=='m'
        res1=(eval(['amp' num2str(order(1))])+eval(['amp' num2str(order(2))]))/2;
        display('Calculating mean of first two peaks...')
    elseif action(1)=='n'
        display('Performing no action with first two peaks...')
    end;


    %Plot if activated
    for PLOT=1:1
    if Fplot==1
        if L>64
            display('Error: Plotting only supports up to 64 epochs')
            return;
        end;
        close all;figure(99);LABLES3={'X0', 'P1','P2','P3'};LABLES=LABLES3(1:inputNum+1);Pos1=FTimeMode;
        %Determine size of plot
        size1=sqrt(L);rest=size1-fix(size1);
        if rest<=0.5
            sizeA=ceil(size1);
            sizeB=floor(size1);
        elseif rest>=0.5
            sizeA=ceil(size1);
            sizeB=ceil(size1);
        elseif rest==0
            sizeA=size1;
            sizeB=size1;
        end;
        if ~isnan(res2)
            plot_results=res2;
        else
            plot_results=res1;
        end;

        for i=1:L
            hand(i)=sbplot(sizeA,sizeB,i);
            plot(EEG.data(electrode,:,i),'g','LineWidth', 1.0);
            switch inputNum %Fit XTicks to number of input Peaks
                case 1
                    [Position, index]=sort([Pos1 time1TP(i)]);
                case 2
                    [Position, index]=sort([Pos1 time1TP(i) time2TP(i)]);
                case 3
                    [Position, index]=sort([Pos1 time1TP(i) time2TP(i) time3TP(i)]);
            end;       
            for ii=1:length(Position)-1 %Case of double Entrys, would cause crash
                if Position(ii)==Position(ii+1)
                    Position(ii+1)=Position(ii+1)+1;
                end;
            end;

            for ii=1:length(index)
                LABLES2{ii}=LABLES{index(ii)};
            end;
            set(gca,'XTick', Position,'XTickLabel',LABLES2,'XGrid','on', 'YGrid','off','YDir','reverse');
            if ~isnan(plot_results)
                set(get(gca,'title'),'String',['VP' num2str(i) ' Result= ' num2str(plot_results(i))]);
            else
                set(get(gca,'title'),'String',['VP' num2str(i) ' Peak 1= ' num2str(amp1(i))]);
            end;

            %Check if equal YLimits are requested, calculate and set them
            if FEqualYLim==1
                pos_lim=ceil(max(max(EEG.data(electrode,:,:))));
                neg_lim=floor(min(min(EEG.data(electrode,:,:))));
                set(gca,'YLim',[neg_lim pos_lim]);
            end;

            set(99,'Position', [1 -17 1440 801])    
        end;
    end;
    end;

    %Construct Variables and return
    switch inputNum
        case 1
            ALL(vpl).all_times={time1ms};ALL(vpl).all_amps={amp1};
        case 2
            ALL(vpl).all_times={time1ms; time2ms};ALL(vpl).all_amps={amp1; amp2};
        case 3
            ALL(vpl).all_times={time1ms; time2ms; time3ms};ALL(vpl).all_amps={amp1; amp2; amp3};
    end;
    %If Calculation was performed, otherwise return NaN (standard setting)
    if ~isnan(res2)
        ALL(vpl).all_results=res2;
    else
        ALL(vpl).all_results=res1;
    end;
    ALL(vpl).EpochSTD=Abweichung;
    
    %Return Trialnumbers if set
    if FRetTrnr==1
        for EPL=1:size(EEG.epoch,2)
            if iscell(EEG.epoch(EPL).eventtrnr)
                TRNR(EPL)=EEG.epoch(EPL).eventtrnr{1};
            else
                TRNR(EPL)=EEG.epoch(EPL).eventtrnr(1);
            end;
        end;
                
        ALL(vpl).Trnr=TRNR;
    end;
    
end;
return;