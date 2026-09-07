function [ output_args ] = AGF_start_ica( inname, outname, inlocation, outlocation, chanlocation, method)
%Simple function that re-references to common average, looks up chanlocs (input: chlocation as str),
%performs 3 types of ICA if wanted (short, extended, AMICA (does currently not work and needs additional files).


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
   
% RE-REFERENCE
EEG = pop_loadset( 'filename', inname, 'filepath',inlocation);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );
EEG = pop_reref( EEG, [], 'refstate',0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'overwrite', 'on');  


% CHANNELINFO
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'overwrite', 'on'); 
EEG=pop_chanedit(EEG,  'lookup', chanlocation)
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'overwrite', 'on');  


% ICA
if ~strcmp(method, 'amica')
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    if strcmp(method, 'extended')
        EEG = pop_runica(EEG,  'icatype', 'runica', 'dataset',1, 'options',{ 'extended',1});
    else
        EEG = pop_runica(EEG,  'icatype', 'runica', 'dataset',1, 'options',{ 'extended',0});
    end;
else
    [EEG.icaweights, EEG.icasphere, mods ] = runamica12( EEG.data(:,:), 'indir', '/Users/Erdnase/Documents/AlexEEGTest/data/amica/' );%
    EEG = eeg_checkset(EEG); eeglab redraw
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_saveset( EEG,  'filename', outname, 'filepath', outlocation);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);