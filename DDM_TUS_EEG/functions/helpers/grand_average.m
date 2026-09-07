function [ ] = grand_average(inlocation, outlocation, VP, condition, varargin )
%Produces GAs for input VP and Conditions. 
%AGF 2010
%defaults
Foutname=condition;            %Default: outnames = condition

%Check n-args in
nargs = nargin-6;
if nargs > 1 
  if ~(round(nargs/2) == nargs/2)
    error('Odd number of input arguments??')
    return
  end;
end;
%Read arg in
for i = 1:2:length(varargin)
    Param = varargin{i};
    if ~isstr(Param)
      error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
     case 'outname'
         Foutname=varargin{i+1};
    end;
end;

if length(Foutname) ~= length(condition)
    disp('Error: Number of outnames does not match condition!')
    return
end;

for it=1:length(condition)   
    
    for ii=1:length(VP)
        if ii==1
            av={['pb' VP{ii} condition{it} '.set']};
        else
            av=[av {['pb' VP{ii} condition{it} '.set']}];
        end;
    end;
    
    [EEG,com] = pop_grandaverage( av, 'pathname', inlocation); 
    
    EEG = eeg_checkset( EEG );    
    EEG = pop_saveset( EEG,  'filename', ['GA_' Foutname{it} '.set'], 'filepath', outlocation);

    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset( 'filename', ['GA_' Foutname{it} '.set'], 'filepath', outlocation);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = eeg_checkset( EEG );
    EEG = pop_editset(EEG, 'setname', ['GA_' Foutname{it}]);
    %EEG=pop_chanedit(EEG,  'lookup', '/Applications/eeglab7_2_9_20b/plugins/dipfit2.2/standard_BEM/elec/standard_1005.elc', 'eval', 'chantmp = pop_chancenter( chantmp, [],[]);', 'load',{ '/daten/cogneuro/eeglab/braincap_eog.ced', 'filetype', 'autodetect'}, 'nosedir', '+Y');
    load('/Volumes/MobileGoose/SMAC/1auswertung/5pessi/eeg/skript/functions/CHANEL.mat')
    EEG.chanlocs=CHANEL;
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG,  'filename', ['GA_' Foutname{it} '.set'], 'filepath', outlocation);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
end
