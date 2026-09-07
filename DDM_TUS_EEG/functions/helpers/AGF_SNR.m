function [SNR] = AGF_SNR(inname, inpath, varargin)
%Plots figures that compare two sessions in one subject at the same electrode.

%Input Arguments
% inname        =   has to be an array of inpput names of filenames (multiple of two if two separate files are used)
% inlocation    =   folder for input files (all have to be in one folder)
% channel       =   channel lable to use
% SNR           =   two timewindows which will be used as signal (first) to
%                   noise (second) calculation (example: {30:50; 90:100})
% Filter        =   Filters the data before measuring SNR. Has to be
%                   entered in the following way: {'bp'; [42; 50]}. First arguments can be: 
%                   'notch', LP (low pass), HP (high pass) or BP (bandpass, required two frequencies). 

%Output Arguments
%SNR            =   Signal to noise ration in the requestet timewindow for all datasets seperately

%SET DEFAULT VALUES
nargs = nargin-2;
Fchannel='FCz';     %Default channel is FCz
FSNR={10:110; -150:-50};
Ffilt=[];           %Default = no filter

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
    case 'channel'
         Fchannel=varargin{i+1}
    case 'snr'
         FSNR=varargin{i+1};
     case 'filter'
         Ffilt=varargin{i+1};
         if (strcmp(Ffilt{1},{'lp'}) || strcmp(Ffilt{1},{'hp'}) || strcmp(Ffilt{1},{'notch'})) && length(Ffilt{2})>1
             disp('Error: For low and high pass filters only one argument is allowed');return
         elseif (strcmp(Ffilt{1},{'lp'}) || strcmp(Ffilt{1},{'hp'}) || strcmp(Ffilt{1},{'notch'})) && length(Ffilt{2})==1
             disp(['Performing ' Ffilt{1} '-filter at frequency ' num2str(Ffilt{2}) 'Hz before SNR calculation.']);
         elseif strcmp(Ffilt{1},{'bp'}) && length(Ffilt{2})~=2
             disp('Error: For band pass filters two arguments are required');return
         elseif strcmp(Ffilt{1},{'bp'}) && length(Ffilt{2})==2
             disp(['Performing ' Ffilt{1} '-filter from ' num2str(Ffilt{2}(1)) ' to ' num2str(Ffilt{2}(2)) 'Hz before SNR calculation.']);
         else
             disp('Error: Unknown Filter settings.');return
         end;
    end;
end; %Get Variables from varargin array

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for c=1:length(inname)
    %LOAD DATA and get curves for plotting
    EEG = pop_loadset( 'filename', inname{c}, 'filepath', inpath);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = eeg_checkset( EEG );
    epochs(c)=EEG.trials;
    %Search channel by name
    CH=find(strcmp({EEG.chanlocs.labels},Fchannel));
    
    %if filter is set: Filter EEG before calculations
    if ~isempty(Ffilt)
        switch Ffilt{1}
            case 'lp'
                EEG = pop_eegfilt( EEG, 0, Ffilt{2}, [], [0], 0, 0, 'fir1', 0);EEG = eeg_checkset( EEG );
            case 'hp'
                EEG = pop_eegfilt( EEG, Ffilt{2}, 0, [], [0], 0, 0, 'fir1', 0);EEG = eeg_checkset( EEG );
            case 'bp'
                EEG = pop_eegfilt( EEG, Ffilt{2}(1), Ffilt{2}(2), [], [0], 0, 0, 'fir1', 0);EEG = eeg_checkset( EEG );
            case 'notch'
                disp('Notch filtering leads to errors in current matlab, not supported ATM'); return
        end;
    end;
    
    TW=[min(EEG.times):max(EEG.times)];
    dat(c,:)=mean(EEG.data(CH,:,:),3); %Average over epochs

    
    if ~isempty(FSNR)
        T_S=[find(EEG.times==min(FSNR{1})):find(EEG.times==max(FSNR{1}))]; %timewindow for signal
        T_N=[find(EEG.times==min(FSNR{2})):find(EEG.times==max(FSNR{2}))]; %timewindow for noise
        
        signal_std{c}=squeeze(std(EEG.data(CH,T_S,:)));
        noise_std{c}=squeeze(std(EEG.data(CH,T_N,:)));
        
        signal_rms{c}=squeeze(sqrt(mean(EEG.data(CH,T_S,:).^2)));
        noise_rms{c}=squeeze(sqrt(mean(EEG.data(CH,T_N,:).^2)));
        
        SN_std(c)=mean(squeeze(std(EEG.data(CH,T_S,:)))./squeeze(std(EEG.data(CH,T_N,:)))); %ratio as average of std of signal / std of noise in each epoch
        SN_rms(c)=mean(squeeze(sqrt(mean(EEG.data(CH,T_S,:).^2)))./squeeze(sqrt(mean(EEG.data(CH,T_N,:).^2)))); %ratio as average of RMS of signal / RMS of noise in each epoch
    end;
    ALLEEG = pop_delset( ALLEEG, [1] );

end; 

SNR.SN_std=SN_std;
SNR.SN_rms=SN_rms;