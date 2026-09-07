%function [ ] = AGF_contrasts( loadpath )
%Function that calculates contrasts and separates datasets from regression
%output of SMAC study (i.e. Genotype x Medication). 
%Requires 'AGF_epoch_by_number', 'easy_topo' and 'differencewave' functions to be present.
%AGF, 2013
%__________________________________________________________________________
%Input:
%'loadpath'     -   Folder with datasets to use. Will automatically use all
%                   available datasets and create subfolders with grafix.
%__________________________________________________________________________
%Index groups by gene and medication
loadpath='/Volumes/AGF work/SMAC/1auswertung/5pessi/eeg/model/data/GLM201/';

SP=1:15; SV=16:30; LP=31:46; LV=47:62;
filez=dir([loadpath '*.set']);

%Check if all subfoldes exist and create them if neccessary
folds={'data group';'ANOVA gene-medication';'contrast gene';'contrast medication'};
for c = 1 : length(folds)
    if exist([loadpath folds{c}])==7
        disp(['Folder: "' folds{c} '" is already present. Data might be overwritten if names collide.'])
        if exist([loadpath folds{c} '/grafix'])~=7
            mkdir([loadpath folds{c}], 'grafix')
        end;
    else
        disp(['Folder: "' folds{c} '" is not present - will be created'])
        mkdir([loadpath folds{c}], 'grafix')
        mkdir([loadpath folds{c}], 't-map')
    end;
end;

%Step 1: Create separate datasets for genetic x medication and all saline and all medication
indi={SP; SV; LP; LV; [SP LP]; [SV LV]};
names={'SS saline';'SS citalopram';'LL saline';'LL citalopram';'All saline';'All citalopram'};
brzk
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for c = 1 : length(filez)
    %Load first dataset that contains all groups
    EEG = pop_loadset( 'filename', filez(c).name, 'filepath', loadpath);
    EEG = AGF_epoch_by_number(EEG, 1:62);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
    
    for cc = 1 : length(indi)
        %Extract all groups separately and save this
        EEG = AGF_epoch_by_number(EEG, indi{cc});
        savename=[names{cc} ' ' EEG.filename]; %Remember dataset name
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off', 'setname', [names{cc} EEG.setname(22:end)]);
        EEG = pop_saveset( EEG, 'filename', savename, 'filepath', [loadpath folds{1}]);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        ALLEEG = pop_delset( ALLEEG, [2] );
        [EEG ALLEEG CURRENTSET] = eeg_retrieve(ALLEEG,1);
    end;
    ALLEEG = pop_delset( ALLEEG, [1] );
end;
%Plot results
easy_topo({[loadpath folds{1} '/']}, [loadpath folds{1} '/grafix/'], 50)
%%

%Step 2: Create contrast for medication within and over groups
%retrieve indices of all created files
iSP=dir([loadpath folds{1} '/*SS saline*.set']);
iSV=dir([loadpath folds{1} '/*SS cit*.set']);
iLP=dir([loadpath folds{1} '/*LL saline*.set']);
iLV=dir([loadpath folds{1} '/*LL cit*.set']);
iAP=dir([loadpath folds{1} '/*All saline*.set']);
iAV=dir([loadpath folds{1} '/*All cit*.set']);

contrastnames1={'contrast t-map all medication', 'contrast t-map SS medication', 'contrast t-map LL medication'};
contrastnames2={'contrast corrected all medication', 'contrast corrected SS medication', 'contrast corrected LL medication'};
indinames={iAP;iAV;iSP;iSV;iLP;iLV};

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for cc = 1 : 3
    for c = 1 : length(iSP)
        %Load first dataset that contains all groups
        EEG = pop_loadset( 'filename', indinames{cc}(c).name, 'filepath', [loadpath folds{1}]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
        EEG = pop_loadset( 'filename', indinames{cc+1}(c).name, 'filepath', [loadpath folds{1}]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 2 );
        
        [cor_d, diff_d, t_d] = pop_differencewaves(ALLEEG(1), ALLEEG(2), 'method','p');
        
        %Save corrected differencewave
        if size(cor_d,3)==2
            EEG = pop_selectevent( EEG,  'epoch', [1:2], 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
        end;
        savename2=[contrastnames2{cc} ' ' filez(c).name]; %Remember dataset name
        EEG.data=cor_d;
        EEG = eeg_checkset( EEG );
        EEG = pop_editset(EEG, 'setname', savename2(1:end-4));
        EEG = pop_saveset( EEG, 'filename', savename2, 'filepath', [loadpath folds{4} '/']);
        
        %Save t-statistic as EEG Data
        if size(t_d,3)==2 && size(cor_d,3)~=2
            EEG = pop_selectevent( EEG,  'epoch', [1:2], 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
        end;
        savename1=[contrastnames1{cc} ' ' filez(c).name]; %Remember dataset name
        EEG.data=t_d;
        EEG = eeg_checkset( EEG );
        EEG = pop_editset(EEG, 'setname', savename2(1:end-4));
        EEG = pop_saveset( EEG, 'filename', savename2, 'filepath', [loadpath folds{4} '/t-map/']);
        
        ALLEEG = pop_delset( ALLEEG, [1 2] );
    end;
end;


%%
%Step 3: Create contrast for genetics
%retrieve indices of all created files
iSP=dir([loadpath folds{1} '/*SS saline*.set']);
iSV=dir([loadpath folds{1} '/*SS cit*.set']);
iLP=dir([loadpath folds{1} '/*LL saline*.set']);
iLV=dir([loadpath folds{1} '/*LL cit*.set']);

contrastnames1={'contrast t-map genetics in saline', 'contrast t-map genetics in medication'};
contrastnames2={'contrast corrected genetics in saline', 'contrast corrected genetics in medication'};
indinames={iSP;iLP;iSV;iLV};

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for cc = 1 : 2
    for c = 1 : length(iSP)
        %Load first dataset that contains all groups
        EEG = pop_loadset( 'filename', indinames{cc}(c).name, 'filepath', [loadpath folds{1}]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
        EEG = pop_loadset( 'filename', indinames{cc+1}(c).name, 'filepath', [loadpath folds{1}]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 2 );
        
        [cor_d, diff_d, t_d] = pop_differencewaves(ALLEEG(1), ALLEEG(2), 'method','p');
        
        %Save corrected differencewave
        if size(cor_d,3)==2
            EEG = pop_selectevent( EEG,  'epoch', [1:2], 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
        end;
        savename2=[contrastnames2{cc} ' ' filez(c).name]; %Remember dataset name
        EEG.data=cor_d;
        EEG = eeg_checkset( EEG );
        EEG = pop_editset(EEG, 'setname', savename2(1:end-4));
        EEG = pop_saveset( EEG, 'filename', savename2, 'filepath', [loadpath folds{3} '/']);
        
        %Save t-statistic as EEG Data
        if size(t_d,3)==2 && size(cor_d,3)~=2
            EEG = pop_selectevent( EEG,  'epoch', [1:2], 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
        end;
        savename1=[contrastnames1{cc} ' ' filez(c).name]; %Remember dataset name
        EEG.data=t_d;
        EEG = eeg_checkset( EEG );
        EEG = pop_editset(EEG, 'setname', savename2(1:end-4));
        EEG = pop_saveset( EEG, 'filename', savename2, 'filepath', [loadpath folds{3} '/t-map/']);
        
        ALLEEG = pop_delset( ALLEEG, [1 2] );
    end;
end;
