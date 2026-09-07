function [  ] = differencewaves( inname, inname2, outname, inlocation, inlocation2, outlocation, varargin)
%function that creates differencewaves of the averaged ERPs of two input
%datasets. Inputs must be before averaging! Function can automaticly
%equalize the epochnumbers if wanted (random trials are picked)... Function will always
%subtract dataset2 from dataset1 (i.e. if 1 is more negative, result stays negative).
%
%
%Optional inputs
%
%'setname' = any string      | Name for the new dataset.
%'modus' = 'equal'           | Equalizes epochnumbers of both input datasets to the smaller one's number.
%'baserem' = [min_ms max_ms] | Baseline latency range in milliseconds. [] for whole epoch (demean). 
%'plot' = 0/channumber       | Creates a simple plot of the differencewaveform at the channel specified after plot and pauses.
%'average' = 0/1             | Datasets are averaged at the end (default is 0 = no averaging over epochs).
%'method' = ['subtract'/'p'] | Defines the method of creating the difference-datasets: 
%                            | 'subtract' just calculates the difference between both datasets.
%                            | 'p' calculates a paired ttest at significancelevel 'alpha'. Nonsignificant datapoints are set to 0.
%                            | 'fdr' uses FDR correctionand resets nonsignificant datapoints to zero (alpha can also be set).
%'fdr' = ['pdep' or 'dep']   | 'pdep' original B&H method, 'dep' B&Y method (valid for all data) is used.
%'alpha' = 0 - 1             | Significance level for tests (default: 0.05)
%'inv' = [1 1]               | Factor to multiply data with: Example [1 -1] inverts dataset 2.
%
%
%AGF, 2010

%defaults
Fmodus='';          %Default: Does not equalize Epochnumbers.
Fsetname=outname;   %Default: Setname is the same as the outname.
Fbaserem=NaN;       %Default: No baseline correction.
Fplot=0;            %Default: No simple plot at the end.
Faverage=0;         %Default: Dataset will not be averaged.
Falpha=0.05;        %Default: Significance level = 0.05.
Fmethod='subtract'; %Default: Datasets be subtracted.
Finv=[1 1];         %Default: No Dataset is inverted.
Ffdr='pdep';        %Default FDR method valid for positive dependency. 


%Check n-args in
nargs = nargin-6;
if nargs > 1 
  if ~(round(nargs/2) == nargs/2)
    error('Odd number of input arguments??')
    return
  end;
end;
%Read arg in
for i = 1:2:length(varargin) %Get input parameters
    Param = varargin{i};
    if ~isstr(Param)
      error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
     case 'modus'
         Fmodus=varargin{i+1};
     case 'baserem'
         Fbaserem=varargin{i+1};
     case 'setname'
         Fsetname=varargin{i+1};   
     case 'plot'
         Fplot=varargin{i+1};
     case 'average'
         Faverage=varargin{i+1};
     case 'alpha'
         Falpha=varargin{i+1};
     case 'method'
         Fmethod=varargin{i+1};
     case 'fdr'
         Ffdr=varargin{i+1};
         if strcmp(Ffdr, 'dep')
             disp('FDR correction method is set to "dep" using Benjamini & Yekutieli method, valid for all dependency structures.')
         end;
     case 'inv'
         Finv=varargin{i+1};
         if size(Finv,2)~=2
             disp('Error: If inv is set, it has to contain one factor for dataset 1 and another one for dataset 2!')
             return
         end;
         disp(['Multiplying dataset 1 with ' num2str(Finv(1)) ' and dataset 2 with ' num2str(Finv(2)) ' !'])
    end;
end;

%Load datasets and get epochnumbers
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset( 'filename', [inname '.set'], 'filepath', inlocation);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = pop_loadset( 'filename', [inname2 '.set'], 'filepath', inlocation2);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EP1=size(ALLEEG(1).epoch,2);EP2=size(ALLEEG(2).epoch,2);
ALLEEG(1).data=ALLEEG(1).data*Finv(1);
ALLEEG(2).data=ALLEEG(2).data*Finv(2);


%If modus is set, Equalize Epochnumbers
if strcmp(Fmodus,'equal')
    if EP1>EP2
        disp(['Reducing Epochnumber to ' num2str(EP2) ' random Epochs as in Dataset ' inname2])
        Index=randperm(EP1);Index=Index(1:EP2);
        EEG = pop_selectevent( ALLEEG(1),  'epoch', Index, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'overwrite', 'on', 'gui', 'off');
    elseif EP2>EP1
        disp(['Reducing Epochnumber to ' num2str(EP1) ' random Epochs as in Dataset ' inname])
        Index=randperm(EP2);Index=Index(1:EP1);
        EEG = pop_selectevent( ALLEEG(2),  'epoch', Index, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'overwrite', 'on', 'gui','off');
    else
        disp(['Both Dataests have the same Epochnumber of ' num2str(EP1)])
    end;
else
    if EP1==EP2
        disp('*****Epochnumbers are identical!******')
    else
        disp('*****Modus not set to "equal" and datasets do not have the same epochnumber!******')
        return
    end;
end;

%If baserem is set, perform baserem before creating Differencewaves
if isempty(Fbaserem)
    disp('De-meaning the whole epoch.')
    ALLEEG(1) = pop_rmbase( ALLEEG(1), [],1:size(ALLEEG(1).data,2));
    ALLEEG(2) = pop_rmbase( ALLEEG(2), [],1:size(ALLEEG(1).data,2));
elseif ~isnan(Fbaserem)
    disp(['Removing Baseline from ' num2str(Fbaserem(1)) 'ms to ' num2str(Fbaserem(2)) 'ms.'])
    ALLEEG(1) = pop_rmbase( ALLEEG(1), Fbaserem);
    ALLEEG(2) = pop_rmbase( ALLEEG(2), Fbaserem);
end;

%Average Datasets
if Faverage==1
    EEG = pop_grandaverage(ALLEEG(1), 'datasets', 1);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'overwrite', 'on', 'gui', 'off','setname',Fsetname); 
    EEG = pop_grandaverage(ALLEEG(2), 'datasets', 1);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'overwrite', 'on', 'gui', 'off','setname',Fsetname);
end;

%Generate Differencewaveform and save new dataset
%Display choosen parameters
if strcmp(Fmethod,'subtract')
    disp('Dataset1 will be subtracted from dataset2...')
    %Subtract datasets
    data=ALLEEG(1).data-ALLEEG(2).data;

elseif strcmp(Fmethod,'p')
    disp(['Performing paired ttest of dataset1 and dataset2 at alpha: ' num2str(Falpha) ' all other datapoints will be set to zero.'])
    data=ALLEEG(1).data-ALLEEG(2).data;
    %Calculat n-ttests.
    for c1=1:size(ALLEEG(1).data,1)
        for c2=1:size(ALLEEG(1).data,2)
            [h(c1,c2),p(c1,c2)]=ttest(ALLEEG(1).data(c1,c2,:), ALLEEG(2).data(c1,c2,:), Falpha);
        end;
    end;
    %Reject data at nonsignificant points.
    for c1=1:size(ALLEEG(1).data,1)
        for c2=1:size(ALLEEG(1).data,2)
            if ~h(c1,c2)
                data(c1,c2,:)=0;
            end;
        end;
    end;
    
elseif strcmp(Fmethod,'fdr')
    disp(['Performing FDR correction of dataset1 and dataset2 at alpha: ' num2str(Falpha) ' all other datapoints will be set to zero.']) 
    data=ALLEEG(1).data-ALLEEG(2).data;
    %Calculat n-ttests.
    for c1=1:size(ALLEEG(1).data,1)
        for c2=1:size(ALLEEG(1).data,2)
            [h(c1,c2),p(c1,c2)]=ttest(ALLEEG(1).data(c1,c2,:), ALLEEG(2).data(c1,c2,:), 0.05);
        end;
    end;    
    %CALL FDR
    [p_adj]=fdr_bh(p, Falpha, 'no', Ffdr);
    %Reject data at nonsignificant points.
    for c1=1:size(ALLEEG(1).data,1)
        for c2=1:size(ALLEEG(1).data,2)
            if ~p_adj(c1,c2)
                data(c1,c2,:)=0;
            end;
        end;
    end;
else
    disp(['Input argument `method´ set to invalid value: ' Fmethod '.'])
    return
end;

%Save
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'retrieve',1, 'study',0); 
EEG=eeg_checkset(EEG);
EEG = pop_editset( EEG, 'setname',  outname );
EEG.data=data;
EEG = pop_saveset( EEG,  'filename', [outname '.set'], 'filepath', outlocation);

%Optional plotting for inspection
if Fplot~=0
    plot(mean(squeeze(mean(ALLEEG(1).data(Fplot,:,:),3)),1));
    legend(Fsetname);
    set(gca,'YGrid','off','YDir','reverse');
    pause
end;