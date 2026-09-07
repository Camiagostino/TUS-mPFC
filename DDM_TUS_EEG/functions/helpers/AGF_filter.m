function [data] = AGF_filter( inname, outname, inlocation, outlocation, low, high, notch, setname )
% SMAC Filter for EEG import
% load data into eeglab
if inname(end-3:end)=='vhdr'
    inlocation
    inname
    EEG = pop_loadbv(inlocation, inname);
else
    EEG = pop_loadset( 'filename', inname, 'filepath', inlocation);
end;    
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'overwrite', 'on');  
EEG = pop_editset(EEG, 'setname', setname);
EEG = eeg_checkset( EEG );

% HighPass filter if set      
if ~isempty(high)
    EEG = pop_eegfilt( EEG, high, 0, [], [0], 0, 0, 'fir1', 0);
    %EEG = pop_eegfilt( EEG, high, 0, [], [0], 'firtype', 'fir1');
    %EEG = pop_iirfilt( EEG, high, 0, [], [1]);
    EEG = eeg_checkset( EEG );
end;

% LowPass filter if set
if ~isempty(low)
    EEG = pop_eegfilt( EEG, 0, low, [], [0], 0, 0, 'fir1', 0);
    %EEG = pop_eegfilt( EEG, 0, low, [], [0], 'firtype', 'fir1');
    %EEG = pop_iirfilt( EEG, 0, low, [], [1]);
    EEG = eeg_checkset( EEG );
end;

%Notch filter if set
if ~isempty(notch)
    EEG = pop_iirfilt( EEG, notch(1), notch(2), [], [1]);
    EEG = eeg_checkset( EEG );
end;

%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'overwrite', 'on');  
EEG = eeg_checkset( EEG );

if ~isempty(outlocation)
    EEG = pop_saveset( EEG,  'filename', outname, 'filepath', outlocation);
end;

data=EEG;
return;

