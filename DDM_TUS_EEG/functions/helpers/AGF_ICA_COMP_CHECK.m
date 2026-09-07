function [ outp, add_info, Odata ] = AGF_ICA_COMP_CHECK( indata, trial_ind, time_w, varargin)
%%
%Function that does something with ICA components.
%
%INPUT:
%'indata'       either EEG structure or path to eeg file that will be
%               openend.
%'trial_ind'    indexes the trials of interest.
%'time_w'       time window (ms) that should be averaged over.
%
%OPTIONAL INPUT:
%
%'TF'           Activates regression in frequency domain (set to string:
%               cell array {4 8} for theta power). Default: {};
%'disp'         Toggle display of percentage [0/1].
%'aleph'        alpha value for the analysis (default: .05)
%'adj'          Adjust p values. Values 'bonf' | 'bonf2' for Bonferroni
%               correction over electrode | time and electrode.
%'bin_size'     Activates multiple options in which the component is 
%               tested. This size here defines the window length.
%'stepsize'     Performs the test in every n-th datapoint (default =
%               length(time_w) -> only one test.
%'return_act'   [0/1] Will return activity for conditions trial_ind (only works
%               if index has two unqiue numbers).
%'increaseR2'   Will first select components and then remove all components
%               that do not increase the backprojection R2 of the best component by at least this factor (1=any
%               increase, 1.05 at least 5%, 2 = at least double...)
%'otherdomain'  [0-1] will also remove components that in the respective
%               other domain (time when frequency is selected) are not significant at the
%               input level. Set = 1 to use same p-value threshold as initial selection,
%               every other number between 0 and 1 will use this p-value.
%'spec_elec'    Numbers of electrodes that should be outliers (all other
%               components are excluded).
%'spec_alpha'   Alpha value for Grubbs Test that determines if the correct
%               electrodes are outliers.
%'include'      Component numbers that are always included.
%'exclude'      Component numbers that are always excluded.
%'n_std'        Number of STD that electrodes should be away (=number). 
%'test_window'  Can be set to a (smaller) value than the action time_w, therefore longer plots
%               are generated and criteria for comp selection are only evaluated in a subset. 
%'standardize'  [0/1] if 1, ICA activity will be zscored per component in the overall timewindow.
%               This will also lead to standardized Beta weights, comparable accross subjects.
%'frequency'    Can be cell array that includes any two numbers {3.5 40}. Will return
%               activity and spectral plot. 
%'Freq_Compress'Will return activity collapsed over these frequencies...
%'freq_base'    Baseline for TF transformation (in ms of original data)
%'TF_stat_time' Time_window to perform between condition statistics for
%               frequencies (ms).
%'TF_window'    Time Window to return TF activity in (default -1000 ms to 1000 ms)
%'BP_electrode' Function will return backprojection at electrode of highest
%               activity and the electrode specified here (can be an array, i.e. [13 22 30].
%'perc_base'    Translates activity into % signal change from a baseline
%               (akin to TF analyses).
%'base'         Specify certain periods per single trial as baseline (has
%               to be an 2D array of length (trial_ind) where column one =
%               first timepoint and column 2 = second timepoint.
%'model n'      Provide a model to perform robust regression in. Model has
%               to have the following structure = {{trialindex} {predictors} {Y} {name}}
%               trialindex = trials to be included relative to input 'trial_ind'.
%               {name} has to be a string.
%               Any number of models can be provided via input 'model 1' 'model 2' 'model 3' etc.
%
%RETURNS:
%'outp'         Contains in the following order: 
%                   1. Number if significant component
%                   2. F value
%                   3. R2 value
%                   4. p value
%                   5. Time of best fit (if bins are activated,  otherwise NaN).
%'add_info'     Contains the input EEG structure without data (can be saved 
%               for plotting of components etc.) and settings for the
%               analysis.
%'Odata'        Struct array of the activity of the components in the best
%               fitting time window (only when prompted).
%
%Set defaults
%get datapoints from milliseconds
Ftime_w = [find(indata.times==time_w(1)):find(indata.times==time_w(2))];
TF           = '0';             % Default: Work in time domain. 
Fcut_p       = 0.05;            % Default alpha value.
Fadj         = 'bonf';          % Default: Bonferroni for number of components is used.
Freturn_act  = 0;               % Default: Do not return ICA component activity.
Fspec_elec   = NaN;             % Default: No electrode criterion.
Fgrubbsalpha = 0.05;            % Default: Value for Grubbs test is 0.05
Ftest_window = Ftime_w;         % Default: test is the same as time window
Fstandardize = 0;               % Default: Do not zscore ICA activity.
Ffrequency   = {};              % Default: Do not return time-frequency power envelopes.
Ffreq_base   = [-1000 -500];    % Default: No baseline set for TF analysis.
FTF_window   = [-1000 1000];    % Default: Timewindow for TF analysis is -1000 to 1000 ms
FBP_electrode= 30;              % Default: Electrode 30 activity is returned as backprojection.
Fmodel       = {};              % Default: No model used, no single datapoint regression.
Finclude     = [];              % Default: No component is always included
Fexclude     = [];              % Default: No component is awlays excluded.
FincreaseR2  = [];              % Default: Do not select components on R2 BP increase
nmodels      = 0;               % Default: No additional models are specified
Fdisp        = 1;               % Default: Display is activated.
Fotherdomain = NaN;             % Default: Do not assess information in the other domain (frequency vs time) 
Fperc_base   = 0;               % Default: Do not perform transposition into % baseline activity
Fbase        = 0;               % Default: No baseline specified.
Fspeed       = 1;               % Default: Extended more includes full spectrum plot
FStepSizeI   = NaN;
FbinsI       = NaN;
%%
for hide = 1:1 %Get Variables from varargin array
    nargs = nargin-3;
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
            case 'adj'
                Fadj=lower(varargin{i+1});
            case 'perc_base'
                Fperc_base=lower(varargin{i+1});
                if isnumeric(Fperc_base)
                    if Fperc_base ~= 0 & Fperc_base ~= 1
                        disp(['Unrecognized value ' num2str(Fperc_base) ' for baseline input. has to be either 1 (active) or 0.']); return;
                    end;
                else
                    disp(['Unrecognized value ' Fperc_base ' for baseline input. has to be either 1 (active) or 0.']); return;
                end;
            case 'base'
                Fbase=lower(varargin{i+1});
                if size(Fbase,2)~=2
                    disp('Error: Baseline matrix has to be 2D with baselines over trials in columns...'); return;
                end;
            case 'aleph'
                Fcut_p=varargin{i+1};
            case 'bin_size'
                FbinsI=varargin{i+1};
            case 'include'
                Finclude=varargin{i+1};
            case 'exclude'
                Fexclude=varargin{i+1};
            case 'stepsize'
                FStepSizeI=varargin{i+1};
            case 'disp'
                Fdisp=varargin{i+1};
            case 'return_act'
                Freturn_act=varargin{i+1};
            case 'spec_elec'
                Fspec_elec=varargin{i+1};
            case 'increaser2'
                FincreaseR2=varargin{i+1};
            case 'spec_alpha'
                Fgrubbsalpha=varargin{i+1};
            case 'standardize'
                Fstandardize=varargin{i+1};
            case 'bp_electrode'
                FBP_electrode=varargin{i+1};
            case 'otherdomain'
                Fotherdomain=varargin{i+1};
                if Fotherdomain > 1 | Fotherdomain <= 0
                    disp('Error: input otherdomain has to be between 0 and 1 and defines the p value threshold for rejection. Stopping...');return;
                end;
            case 'frequency'
                Ffrequency=varargin{i+1};
                disp(['Returns power envelope for frequency/frequencies: ' Ffrequency{:}])
            case 'freq_base'
                Ffreq_base=varargin{i+1};
            case 'tf_stat_time'
                FTF_stat_time=varargin{i+1};
            case 'tf_window'
                FTF_window=varargin{i+1};
            case 'tf'
                TF=varargin{i+1};
            case ['model ' num2str(nmodels+1)]
                nmodels=nmodels+1;
                Fmodel{nmodels}=varargin{i+1};
                if size(Fmodel{nmodels},2)~=4
                    disp('Models have to have 4 entrys: 1st = trials, 2nd = predictors, 3rd = observations,  4th = name.');return;
                end;
            case 'test_window'
                Ftemp=varargin{i+1};
                Ftest_window=[find(indata.times==Ftemp(1)):find(indata.times==Ftemp(2))];
                disp(['Tests will be performed from DP ' num2str(Ftest_window(1)) ' to ' num2str(Ftest_window(end)) ' and activity returned for DP ' num2str(Ftime_w(1)) ' to ' num2str(Ftime_w(end)) '. Timelocking event at DP: ' num2str(find(indata.times==0)) '.'])
            otherwise
                disp(['Unknown argument ' varargin{i} '...'])
                pause
        end;
    end; 
end;

%%
%Set simple variables and pre allocate arrays
Fbins           = length(Ftest_window);                 % Default: Average over length of time-window
FStepSize       = length(Ftest_window);                 % Default: Average over length of time-window
n_comp          = size(indata.icaweights,1);            % just the number of ICA components
comp_timebest   = NaN;
vas             = unique(trial_ind);                    % Find if regressor is categorical or not
I1              = find([trial_ind]==vas(1));            % Unique entrys assuming binary categorical regressor
I2              = find([trial_ind]==vas(2));
ranks           = nan(n_comp,length(1 : FStepSize : length(Ftest_window)-1),12);  %ranks of the sorted components at the end


if ~isnan(FbinsI)
    Fbins=FbinsI;
end;
if ~isnan(FStepSizeI)
    FStepSize=FStepSizeI;
end;
disp(['Bin size is: ' num2str(Fbins) ' | Stepsize is : ' num2str(FStepSize)])
%Add some initial information to output
add_info.Statistic_Timewindow_ms=[indata.times(Ftest_window(1)) indata.times(Ftest_window(end))];
add_info.Return_Timewindow_ms=time_w;
add_info.Bin_size_ms=Fbins*2;
add_info.Step_size_ms=FStepSize*2;
add_info.Output_Electrodes=FBP_electrode;
add_info.Scale_2_Baseline_Percent=Fperc_base;

for c = 1:length(FBP_electrode)
    OutLable{c}=indata.chanlocs(FBP_electrode(c)).labels;
end;
add_info.Output_Labels=OutLable;  
    
if strcmp(TF,'0')
    add_info.Domain_Selection='Time Domain';
else
    add_info.Domain_Selection=TF;
end;
%%
if Fstandardize
    disp('Zscoring component activations...')
    for z = 1 : n_comp
        indata.icaact(z,:,:)=zscore(indata.icaact(z,:,:));
    end;
end;

%%
if Fbase
    if size(Fbase,1)==1   
        FUbase=repmat(Fbase,size(indata.data,3),1);
    else
        FUbase=Fbase;
    end;
    for trial = 1 : size(indata.data,3)
        if mod(FUbase(trial,1),2)
            b1=find(indata.times==FUbase(trial,1)-1);
        else
            b1=find(indata.times==FUbase(trial,1));
        end;
        if mod(FUbase(trial,2),2)
            b2=find(indata.times==FUbase(trial,2)-1);
        else
            b2=find(indata.times==FUbase(trial,2));
        end;
        %Get single trial baseline
        bsl=repmat(mean(indata.data(:,b1:b2,trial),2),1,size(indata.data,2));
        %subtract by baseline
        Data_base(:,:,trial)=indata.data(:,:,trial)-bsl;
    end;
    %indata.data=Data_base;
    %indata.icaact=[];
    %[indata]=eeg_checkset(indata);
    add_info.Individual_Trial_Baseline=[mean(FUbase(:,1)) mean(FUbase(:,2))]
else
    add_info.Individual_Trial_Baseline=NaN;
end;

%%
%Perform security check: Is enough data at the borders of the time-window for the selected bin size?
if isempty(indata.icaact)
    indata.icaact = eeg_getica(indata);
end;
if Ftime_w(1)-ceil(Fbins/2) <= 0 | Ftime_w(end)+floor(Fbins/2) > size(indata.icaact,2)
    disp('Warning! Timewindow + bin size to average is bigger than input data... stopping.')
    return;
end;

%%
add_info.Specific_Electrodes=Fspec_elec;
add_info.Specific_Electrodes_Criterion=Fgrubbsalpha;
FinclComp       = 1:size(indata.icaweights,1);
if ~isnan(Fspec_elec) %If specified, exclude all components that are not outliers at some certain electrodes
    if Fgrubbsalpha<1
        disp(['Using only components where electrodes '  num2str(Fspec_elec) ' are outliers. Alpha for significance: ' num2str(Fgrubbsalpha)])
        REP=1;
        while REP %loop until at least one component fullfills criterion
            del_ind = []; K = 0;
            del_name={};
            for c = 1 : n_comp
                [~, ~, ~, ~, ~, outlier, outlier_num] = grubbs(indata.icawinv(:,c), {indata.chanlocs.labels}, 1, Fgrubbsalpha);
                if ~isempty(outlier_num)
                    if sum(ismember(Fspec_elec, outlier_num(:,1))) == 0 %outliers do not match input electrodes
                        K = K + 1;
                        del_ind(K) = c;
                    end;
                else
                    K = K + 1;
                    del_ind(K) = c;
                end;
            end;
            if length(del_ind) == length(FinclComp)
                Fgrubbsalpha=Fgrubbsalpha+Fgrubbsalpha*0.1; %Increase alpha threshold until one component is found.
            else
                REP=0;
            end;
        end;
        disp(['Removing ' num2str(length(del_ind)) ' component(s) at alpha ' num2str(Fgrubbsalpha) '. These are: ' num2str(del_ind) '  | ' num2str(length(FinclComp)-length(del_ind)) ' components left.'])
        FinclComp(del_ind)=[];
    else
        %Select components based on defined electrode criteria in relation to SD
        disp(['Using only components where electrodes '  num2str(Fspec_elec) ' are at least ' num2str(Fgrubbsalpha) ' sd more active.'])
        del_ind = []; K = 0;
        for c = 1 : n_comp
            outl = find(abs([indata.icaweights(c,:)]-mean([indata.icaweights(c,:)]))>Fgrubbsalpha*std([indata.icaweights(c,:)]));
            if sum(ismember(Fspec_elec, outl)) == 0 %outliers do not match input electrodes
                K = K + 1;
                del_ind(K) = c;
            end;
        end;
        disp(['Removing ' num2str(length(del_ind)) ' component(s). These are: ' num2str(del_ind) '  | ' num2str(length(FinclComp)-length(del_ind)) ' components left.'])
        add_info.ICs_Removed_by_Spec_Electr=FinclComp(del_ind);
        FinclComp(del_ind)=[];
    end;
end;
add_info.n_ICs_Removed_by_Spec_Electr=length(del_ind);

if ~isempty(Finclude)
    FinclComp=unique([FinclComp Finclude]);
end;
if ~isempty(Fexclude)
    FinclComp=setdiff(FinclComp, Fexclude);
end;
add_info.ICs_Manually_Included=Finclude;
add_info.ICs_Manually_Excluded=Fexclude;
%%
%Adjust significance level
add_info.Initial_alpha=Fcut_p;
if strcmp(Fadj, 'bonf')
    Fcut_p = Fcut_p / length(FinclComp);
    disp(['Bonferroni correction is used. Critical p value is: ' num2str(Fcut_p)])
elseif strcmp(Fadj, 'bonf2')
    Fcut_p = Fcut_p / (length(FinclComp)*length(Ftest_window));
    disp(['Bonferroni correction over time (' num2str(length(Ftest_window)) ') and components (' num2str(length(FinclComp)) ') is used. Critical p value is: ' num2str(Fcut_p)])
end;
add_info.P_caluculation=Fadj;
add_info.Critical_p_value=Fcut_p;
%%
%Decide if component activity should be transformed to time-frequency or not?
if ~strcmp(TF,'0')
    disp(['Using ' TF ' band power for IC selection...'])
    [RegData] = AGF_Hilbert(indata.icaact(FinclComp,:,1:length(trial_ind)), Ffreq_base, {TF}, [indata.xmin*1000 indata.xmax*1000], FTF_window(1), FTF_window(2), indata.srate,'divide');
    RegData = squeeze(RegData);
    StartCorr = Ftest_window(1)-find(indata.times==FTF_window(1));
    MsCorr = find(indata.times==FTF_window(1));
else
    RegData = indata.icaact(FinclComp,:,1:length(trial_ind));
    StartCorr = Ftest_window(1);
    MsCorr= 0;
end;

%%
%Run the regression analysis first in the time-window selected (will be repeated for broader time-range if plot information and test time window are not identical)!
disp(' ');fprintf(1,['Processing subject ' indata.setname(1:6) '\n']);
fprintf(1,['Progress: Percent of analysis done:   ']);
count=0;
for bin = 1 : FStepSize : length(Ftest_window)-1
    count=count+1;
    do=floor((bin/(length(Ftest_window)-1))*100);
    if do < 10
        fprintf(2,'\b%d', do);
    else
        fprintf(2,'\b\b%d', do);
    end;
    
    %Define datapoints: Actual point + bin_size!
    tbin = [bin-ceil((Fbins-1)/2)+StartCorr : bin+floor((Fbins-1)/2)+StartCorr];
    
    %Do this for all components
    for c = 1 : length(FinclComp)
        [t1, t2, ~, ~, STATS] = regress(squeeze(mean(RegData(c, tbin, :),2)), [ones(length(trial_ind),1) trial_ind']);
        ranks(c,count,1)=STATS(2);    %remember F value
        ranks(c,count,2)=STATS(1);    %remember R2 value
        ranks(c,count,3)=STATS(3);    %remember p value
        ranks(c,count,4)=t1(2);       %remember B
        ranks(c,count,5)=t2(2,1);     %remember lower 95 CI
        ranks(c,count,6)=t2(2,2);     %remember upper 95 CI
        ranks(c,count,7)  = mean(squeeze(mean(RegData(c, tbin, :),2)));                           %remember activity all
        ranks(c,count,10) = std(squeeze(mean(RegData(c, tbin, :),2)))/sqrt(length(trial_ind)-1);  %remember std cat all
        if length(vas)==2
            ranks(c,count,8)  = mean(squeeze(mean(RegData(c, tbin, I1),2)));                      %remember activity cat 1
            ranks(c,count,9)  = mean(squeeze(mean(RegData(c, tbin, I2),2)));                      %remember activity cat 2
            ranks(c,count,11) = std(squeeze(mean(RegData(c, tbin, I1),2)))/sqrt(length(I1)-1);    %remember std cat 1
            ranks(c,count,12) = std(squeeze(mean(RegData(c, tbin, I2),2)))/sqrt(length(I2)-1);    %remember std cat 2
        end;
    end;
end;
fprintf(2,'\n')

%Find for each component the timepoint of best fit
for c = 1 : length(FinclComp)
    [comp_F(c), comp_timebest(c)]   = max(ranks(c,:,1)); 
    comp_p(c)                       = ranks(c,comp_timebest(c),3);
    comp_R2(c)                      = ranks(c,comp_timebest(c),2);
end;
%Sort components by decending F values
comp_F          = comp_F';
[comp_F,comp_n] = sort(abs(comp_F),'descend');
comp_R2         = comp_R2(comp_n)';
FinclComp       = FinclComp(comp_n);
comp_p          = comp_p(comp_n)';
%Translate this into milliseconds relative to the time-locking event 0
comp_timebestMS = indata.times([StartCorr+comp_timebest(comp_n)+MsCorr])';
%Find significant components below p-value threshold
i_sign_comp     = find(comp_p < Fcut_p);

%%
%Include manually selected components if they do not pass this threshold
if ~isempty(Finclude)
    for c = 1 : length(Finclude)
        find([FinclComp]==Finclude(c)) > i_sign_comp(end)
        if find([FinclComp]==Finclude(c)) > i_sign_comp(end)
            i_sign_comp(end+1)=length(i_sign_comp)+1; %One more significant component
            FinclComp(i_sign_comp(end))=Finclude(1);
        end;
    end;
end;

add_info.Number_of_significant_ICs=length(i_sign_comp);
add_info.Index_of_significant_ICs=FinclComp(i_sign_comp);
%%
add_info.R2_Rejected_ICs=NaN;
add_info.R2_Criterion=FincreaseR2;
add_info.R2_Additional_ICs_Percent_Diference=NaN;
Delete_Index=[];
if ~isempty(FincreaseR2) & length(i_sign_comp) > 1
    disp(['Only components will be kept that increase R2 of the backprojection at ' indata.chanlocs(FBP_electrode(1)).labels ' by at least ' num2str((FincreaseR2-1)*100) '%.'])
    %See if frequency or time domain are criteria...
    %Step 1: Get best component R2 alone:
    retain = setdiff([1:size(indata.icawinv,2)], FinclComp(1));
    BP_ALL = pop_subcomp( indata, retain, 0); 
    BP_ALL=eeg_checkset(BP_ALL);
    RegDat=BP_ALL.data(:,:,1:length(trial_ind));
    ranksBP = nan(length(Ftime_w),2);
    Istat=find(Ftime_w==Ftest_window(1)):find(Ftime_w==Ftest_window(end));
    AC=FBP_electrode(1);
    if ~strcmp(TF,'0') %Freq domain
        warning off
            [Frequencies]=AGF_Hilbert(RegDat(AC, :, :), Ffreq_base, {TF}, [indata.xmin*1000 indata.xmax*1000], FTF_window(1), FTF_window(2), indata.srate,'divide');
        warning on
        tfrex=find(indata.times==FTF_window(1)):find(indata.times==FTF_window(2));
        interval=find(tfrex==Ftime_w(1)):find(tfrex==Ftime_w(end));                        
        for bin = 1 : length(interval) %Rerun regression in frequency domain
            tbin                    = [interval(bin)-ceil((Fbins-1)/2) : interval(bin)+floor((Fbins-1)/2)];
            [~, ~, ~, ~, STATS]     = regress(squeeze(mean(Frequencies(1,1,tbin,:),3)), [ones(length(trial_ind),1) trial_ind']);
            ranksBP(bin,1)          = STATS(2);    %remember F value
            ranksBP(bin,2)          = STATS(1);    %remember R2 value
        end;
    else %Time Domain
        for bin = 1 : length(Ftime_w)
            %Define datapoints: Actual point + bin_size relative to time_window
            tbin                    = [Ftime_w(bin)-ceil((Fbins-1)/2) : Ftime_w(bin)+floor((Fbins-1)/2)];
            [~, ~, ~, ~, STATS]     = regress(squeeze(mean(RegDat(AC, tbin, :),2)), [ones(length(trial_ind),1) trial_ind']);
            ranksBP(bin,1)          = STATS(2);    %remember F value
            ranksBP(bin,2)          = STATS(1);    %remember R2 value
        end;
    end;
    [~,y]=max(ranksBP(Istat,1));
    y=y+Istat(1)-1;
    Best_IC_R2(1) = ranksBP(y,2);

    %Do this again and add the next components to the BP
    for BP_counter = 2:length(i_sign_comp)
        retain = setdiff([1:size(indata.icawinv,2)], [FinclComp(1) FinclComp(BP_counter)]);
        BP_ALL = pop_subcomp( indata, retain, 0); 
        BP_ALL=eeg_checkset(BP_ALL);
        RegDat=BP_ALL.data(:,:,1:length(trial_ind));
        if ~strcmp(TF,'0') %Freq domain
            warning off
                [Frequencies]=AGF_Hilbert(RegDat(AC, :, :), Ffreq_base, {TF}, [indata.xmin*1000 indata.xmax*1000], FTF_window(1), FTF_window(2), indata.srate,'divide');
            warning on
            tfrex=find(indata.times==FTF_window(1)):find(indata.times==FTF_window(2));
            interval=find(tfrex==Ftime_w(1)):find(tfrex==Ftime_w(end));                        
            for bin = 1 : length(interval) %Rerun regression in frequency domain
                tbin                    = [interval(bin)-ceil((Fbins-1)/2) : interval(bin)+floor((Fbins-1)/2)];
                [~, ~, ~, ~, STATS]     = regress(squeeze(mean(Frequencies(1,1,tbin,:),3)), [ones(length(trial_ind),1) trial_ind']);
                ranksBP(bin,1)          = STATS(2);    %remember F value
                ranksBP(bin,2)          = STATS(1);    %remember R2 value
            end;
        else %Time Domain
            for bin = 1 : length(Ftime_w)
                %Define datapoints: Actual point + bin_size relative to time_window
                tbin                    = [Ftime_w(bin)-ceil((Fbins-1)/2) : Ftime_w(bin)+floor((Fbins-1)/2)];
                [~, ~, ~, ~, STATS]     = regress(squeeze(mean(RegDat(AC, tbin, :),2)), [ones(length(trial_ind),1) trial_ind']);
                ranksBP(bin,1)          = STATS(2);    %remember F value
                ranksBP(bin,2)          = STATS(1);    %remember R2 value
            end;
        end;
        [~,y]=max(ranksBP(Istat,1));
        y=y+Istat(1)-1;
        Best_IC_R2(BP_counter) = ranksBP(y,2);
        Perc_Increase(BP_counter) = (Best_IC_R2(BP_counter)/Best_IC_R2(1)-1)*100;
        if Best_IC_R2(BP_counter)>Best_IC_R2(1)
            if Perc_Increase(BP_counter) < (FincreaseR2-1)*100
                if ~isempty(Finclude)
                    if FinclComp(BP_counter) ~= Finclude
                        disp(['IC ' num2str(FinclComp(BP_counter)) ' increases best R2 by ' num2str(Perc_Increase(BP_counter)) '. Lower than criterion - removing...'])
                        Delete_Index=[Delete_Index BP_counter];
                    else
                        disp(['IC ' num2str(FinclComp(BP_counter)) ' kept due to manual selection.'])
                    end;
               else
                   disp(['IC ' num2str(FinclComp(BP_counter)) ' increases best R2 by ' num2str(Perc_Increase(BP_counter)) '. Lower than criterion - removing...'])
                   Delete_Index=[Delete_Index BP_counter];
               end;
            else
                disp(['IC ' num2str(FinclComp(BP_counter)) ' increases best R2 by ' num2str(Perc_Increase(BP_counter)) '%. Component kept.'])
                Best_IC_R2(1) = Best_IC_R2(BP_counter);
            end;
        else
            if ~isempty(Finclude)
                if FinclComp(BP_counter) ~= Finclude
                    disp(['IC ' num2str(FinclComp(BP_counter)) ' decreases best R2 to ' num2str(Perc_Increase(BP_counter)) '% of its initial value - removing...'])
                    Delete_Index=[Delete_Index BP_counter];
                else
                    disp(['IC ' num2str(FinclComp(BP_counter)) ' kept due to manual selection.'])
                end;
           else
               disp(['IC ' num2str(FinclComp(BP_counter)) ' decreases best R2 to ' num2str(Perc_Increase(BP_counter)) '% of its initial value - removing...'])
               Delete_Index=[Delete_Index BP_counter];
           end;
            
        end;
    end;
    i_sign_comp(Delete_Index)=[];
    add_info.R2_Rejected_ICs=FinclComp(Delete_Index);
    add_info.R2_Additional_ICs_Percent_Diference=Perc_Increase(2:end);
    clear ranksBP
end;
add_info.Number_of_R2_Rejected_ICs=length(Delete_Index);
%%
%Reject components in the respective other domain if set
add_info.Other_Domain_p=Fotherdomain;OD_Index=[];OD_Pval=NaN;
if ~isnan(Fotherdomain)
    if Fotherdomain==1
        p_OD = Fcut_p;
        add_info.Other_Domain_Critical_p=Fcut_p;
    else
        p_OD = Fotherdomain;
    end;
    RegDat=indata.icaact(:,:,1:length(trial_ind));
    for OD_counter = 2:length(i_sign_comp)
        AC=FinclComp(i_sign_comp(OD_counter));  %current component
        ranksOD = nan(length(Ftime_w),2);
        Istat=find(Ftime_w==Ftest_window(1)):find(Ftime_w==Ftest_window(end));
        if ~strcmp(TF,'0') %Freq domain for initial selection
            if OD_counter==2;disp(['Only components will be kept that are also significant in the time-domain at p value: ' num2str(p_OD)]);end;
            warning off
                [Frequencies]=AGF_Hilbert(RegDat(AC, :, :), Ffreq_base, {TF}, [indata.xmin*1000 indata.xmax*1000], FTF_window(1), FTF_window(2), indata.srate,'divide');
            warning on
            tfrex=find(indata.times==FTF_window(1)):find(indata.times==FTF_window(2));
            interval=find(tfrex==Ftime_w(1)):find(tfrex==Ftime_w(end));                        
            for bin = 1 : length(interval) %Rerun regression in frequency domain
                tbin                    = [interval(bin)-ceil((Fbins-1)/2) : interval(bin)+floor((Fbins-1)/2)];
                [~, ~, ~, ~, STATS]     = regress(squeeze(mean(Frequencies(1,1,tbin,:),3)), [ones(length(trial_ind),1) trial_ind']);
                ranksOD(bin,1)          = STATS(2);    %remember F value
                ranksOD(bin,2)          = STATS(1);    %remember R2 value
            end;
        else %Time domain for initial selection
            if OD_counter==2;disp(['Only components will be kept that are also significant in the frequency-domain at p value: ' num2str(p_OD)]);end;
            for bin = 1 : length(Ftime_w)
                %Define datapoints: Actual point + bin_size relative to time_window
                tbin                    = [Ftime_w(bin)-ceil((Fbins-1)/2) : Ftime_w(bin)+floor((Fbins-1)/2)];
                [~, ~, ~, ~, STATS]     = regress(squeeze(mean(RegDat(AC, tbin, :),2)), [ones(length(trial_ind),1) trial_ind']);
                ranksOD(bin,1)          = STATS(2);    %remember F value
                ranksOD(bin,2)          = STATS(3);    %remember R2 value
            end;
        end;
        [~,y]=max(ranksOD(Istat,1));
        y=y+Istat(1)-1;
        OD_Pval(OD_counter)=ranksOD(y,2);
        if ranksOD(y,2) <= p_OD
            disp(['IC ' num2str(FinclComp(OD_counter)) ' has p value ' num2str(ranksOD(y,2)) ' < ' num2str(p_OD) '. Keeping component.'])
        else
            disp(['IC ' num2str(FinclComp(OD_counter)) ' has p value ' num2str(ranksOD(y,2)) ' > ' num2str(p_OD) '. Removing component.'])
            OD_Index=[OD_Index OD_counter];
        end;
    end;
    i_sign_comp(OD_Index)=[];
    add_info.Otherdomain_Rejected_ICs=FinclComp(OD_Index);
    add_info.Otherdomain_Rejected_pvals=OD_Pval(2:end);
end;
add_info.Number_of_Otherdomain_Rejected_ICs=length(OD_Index);

%%
if isempty(i_sign_comp)
    disp('No ICA component shows significant effects! Returnin NaN and stopping here...')
    outp = [NaN NaN NaN NaN NaN];
    Odata = NaN;
    add_info = NaN;
    return;
else
    disp([num2str(length(i_sign_comp)) ' components are significant. These are: ' num2str(FinclComp(i_sign_comp)) ])
    outp = [FinclComp(i_sign_comp)' comp_F(i_sign_comp) comp_R2(i_sign_comp) comp_p(i_sign_comp) comp_timebestMS(i_sign_comp)];
end;



%%
if Freturn_act==1
    L_ind  = [1 4 5 9  10 11 12 18 19 20 21 27 28 29 34 35 36 37 42 43 44 45 51 52 56];  % left  scalp electrodes
    R_ind  = [3 6 7 14 15 16 17 23 24 25 26 31 32 33 38 39 40 41 47 48 49 50 54 55 58];  % right scalp electrodes
    %Use information on best components and extract their activity
    if length(vas) > 2
        disp('Index of trials contains more than 2 unique entrys. Cannot return component activation in categorial form.');
    else
        disp(['Using entrys ' num2str(vas(1)) ' and ' num2str(vas(2)) ' to return categorical components activity.'])
    end;
    for ic = 1 : length(i_sign_comp)+2 %loop over significant components and ERP and backprojection of all components
        disp(' ')
        Odata(ic).setname = indata.setname;
        Odata(ic).temporal_smoothing = Fbins*2;
        %%ADJUST FOR SITUATION: IC = 1 -> ERP, IC = 2 -> Backprojection, IC > 2 -> IC of that number
        if ic==1
            Odata(ic).FieldData = ['ERP Activity @ ' indata.chanlocs(FBP_electrode(1)).labels];
            RegDat=indata.data(:,:,1:length(trial_ind));
            SpectDat=indata;SpectDat.data=SpectDat.data(:,:,1:length(trial_ind));
            NELECTRODES=length(FBP_electrode);
            Istat=find(Ftime_w==Ftest_window(1)):find(Ftime_w==Ftest_window(end));%Get index to restrict timewindow for statistics
        elseif ic==2
            Odata(ic).FieldData = ['BP Activity @ ' indata.chanlocs(FBP_electrode(1)).labels];
            %Calculate backprojection of all components
            retain = setdiff([1:size(indata.icawinv,2)], FinclComp(i_sign_comp));
            BP_ALL = pop_subcomp( indata, retain, 0); 
            BP_ALL=eeg_checkset(BP_ALL);
            RegDat=BP_ALL.data(:,:,1:length(trial_ind));
            SpectDat=BP_ALL;SpectDat.data=SpectDat.data(:,:,1:length(trial_ind));
            NELECTRODES=length(FBP_electrode);
            Istat=find(Ftime_w==Ftest_window(1)):find(Ftime_w==Ftest_window(end));%Get index to restrict timewindow for statistics
        else
            AC=FinclComp(i_sign_comp(ic-2));  %current component
            Odata(ic).FieldData = ['IC No ' num2str(AC)];
            CurrLable='Field1_';
            RegDat=indata.icaact(:,:,1:length(trial_ind)); NELECTRODES=1;
            SpectDat=indata;SpectDat.data=SpectDat.icaact(:,:,1:length(trial_ind));
        end;
        Odata(ic).component_symmetry_p  = [];            
        Odata(ic).component_symmetry_t  = []; 
        Odata(ic).top_lable             = {};
        Odata(ic).top_value_SD          = [];
        Odata(ic).top_value_abs         = [];
        Odata(ic).outliers              = {};
        Odata(ic).Topography            = [];
        Odata(ic).Freq_X_Axis = [FTF_window(1) : 2 : FTF_window(2)];
                
        %COUNT THROUGH ALL OUTPUT ELECTRODES
        for CurrELEC=1:NELECTRODES
            if ic < 3
                AC=FBP_electrode(CurrELEC);
                CurrLable=['Field' num2str(CurrELEC) '_'];
                disp(['Calculating specific information for ' Odata(ic).FieldData ' and output for: ' add_info.Output_Labels{CurrELEC}])
            else
                disp(['Calculating specific information for ' Odata(ic).FieldData])
            end;
            if length(vas) == 2
                Odata(ic).( 'n_cat1')=length(I1);
                Odata(ic).( 'n_cat2')=length(I2);
            end;

            
            %%%%%%%%%%%%GET DATA FOR TIME DOMAIN%%%%%%%%%%%%%%%%%
            %new ranks but only for selected components
            ranks2 = nan(length(Ftime_w),12);
            if Fperc_base%Scale by baseline
                clear basl_ST_RMS
                if isempty(Fbase)
                    disp('Error: To scale by RMS in a baseline, the baseine has to be specified.');return;
                end;
                for trial = 1 : size(indata.data,3)
                    if mod(FUbase(trial,1),2)
                        b1=find(indata.times==FUbase(trial,1)-1);
                    else
                        b1=find(indata.times==FUbase(trial,1));
                    end;
                    if mod(FUbase(trial,2),2)
                        b2=find(indata.times==FUbase(trial,2)-1);
                    else
                        b2=find(indata.times==FUbase(trial,2));
                    end;
                    %Get single trial baseline
                    basl_ST_RMS(trial,:)=mean(sqrt([RegDat(AC,b1:b2,trial)].^2),2);
                end;
                disp('Scaling output to RMS of baseline...');
                RegDat=RegDat./mean(basl_ST_RMS);
            end;
            %Rerun regression over longer time-window if set
            for bin = 1 : length(Ftime_w)
                %Define datapoints: Actual point + bin_size relative to time_window
                tbin                    = [Ftime_w(bin)-ceil((Fbins-1)/2) : Ftime_w(bin)+floor((Fbins-1)/2)];
                %Do this for all components
                [t1, t2, ~, ~, STATS]   = regress(squeeze(mean(RegDat(AC, tbin, :),2)), [ones(length(trial_ind),1) trial_ind']);
                ranks2(bin,1)           = STATS(2);    %remember F value
                ranks2(bin,2)           = STATS(1);    %remember R2 value
                ranks2(bin,3)           = STATS(3);    %remember p value
                ranks2(bin,4)           = t1(2);       %remember B
                ranks2(bin,5)           = t2(2,1);     %remember lower 95 CI
                ranks2(bin,6)           = t2(2,2);     %remember upper 95 CI
                ranks2(bin,7)           = mean(squeeze(mean(RegDat(AC, tbin, :),2)));                              %remember activity all
                ranks2(bin,10)          = std(squeeze(mean(RegDat(AC, tbin, :),2)))/sqrt(length(trial_ind)-1);     %remember std cat all
                if length(vas)==2
                    ranks2(bin,8)       = mean(squeeze(mean(RegDat(AC, tbin, I1),2)));                              %remember activity cat 1
                    ranks2(bin,9)       = mean(squeeze(mean(RegDat(AC, tbin, I2),2)));                              %remember activity cat 2
                    ranks2(bin,11)      = std(squeeze(mean(RegDat(AC, tbin, I1),2)))/sqrt(length(I1)-1);            %remember std cat 1
                    ranks2(bin,12)      = std(squeeze(mean(RegDat(AC, tbin, I2),2)))/sqrt(length(I2)-1);            %remember std cat 2
                end;
            end;
            Odata(ic).([CurrLable 'timecourse1_all'])    =ranks2(:,7);
            Odata(ic).([CurrLable 'timecourse1_all_se']) =ranks2(:,10);
            if length(vas) == 2
                Odata(ic).([ CurrLable 'timecourse1'])   =ranks2(:,8);
                Odata(ic).([ CurrLable 'timecourse2'])   =ranks2(:,9);
                Odata(ic).([ CurrLable 'timecourse1_se'])=ranks2(:,11);
                Odata(ic).([ CurrLable 'timecourse2_se'])=ranks2(:,12);
            else
                Odata(ic).([ CurrLable 'timecourse1'])   =ranks2(:,8);
                Odata(ic).([ CurrLable 'timecourse2'])   =NaN;
                Odata(ic).([ CurrLable 'timecourse1_se'])=ranks2(:,11);
                Odata(ic).([ CurrLable 'timecourse2_se'])=NaN;
            end;
            %Get stats to compare in the initial timewindow to select components
            [x,y]=max(ranks2(Istat,1));
            y=y+Istat(1)-1;
            %Save best activity
            if isnan(Fbins)
                tbin= Ftest_window+Ftime_w(1);
            else
                tbin = [Ftime_w(y)-ceil((Fbins-1)/2) : Ftime_w(y)+floor((Fbins-1)/2)];
            end;
            Odata(ic).([ CurrLable 'Activitiy_at_best_fit'])    = squeeze(mean(RegDat(AC, tbin, :)))';
            Odata(ic).([ CurrLable 'Value_Cat1'])               = mean(squeeze(mean(RegDat(AC, tbin, I1))));
            Odata(ic).([ CurrLable 'Value_Cat1_se'])            = std(squeeze(mean(RegDat(AC, tbin, I1))))/sqrt(length(I1)-1);
            Odata(ic).([ CurrLable 'Value_Cat2'])               = mean(squeeze(mean(RegDat(AC, tbin, I2))));
            Odata(ic).([ CurrLable 'Value_Cat2_se'])            = std(squeeze(mean(RegDat(AC, tbin, I2))))/sqrt(length(I2)-1);
            Odata(ic).([ CurrLable 'Highest_F'])     = x;
            Odata(ic).([ CurrLable 'Highest_R2'])    = ranks2(y,2);
            Odata(ic).([ CurrLable 'Lowest_p'])      = ranks2(y,3);
            Odata(ic).([ CurrLable 'Critical_p'])    = Fcut_p;
            Odata(ic).([ CurrLable 'timepoint_ms'])  = indata.times([Ftime_w(1)+y])';
            %save correlation coefficients and confidence intervals
            Odata(ic).([ CurrLable 'beta_timecourse'])      = ranks2(:,4);
            Odata(ic).([ CurrLable 'beta_timecourse_UCI'])  = ranks2(:,5);
            Odata(ic).([ CurrLable 'beta_timecourse_LCI'])  = ranks2(:,6);
            if ic==1
                TOP=squeeze(mean(mean(RegDat(:, tbin, I1),3),2))'-squeeze(mean(mean(RegDat(:, tbin, I2),3),2))';        %For ERP topo is topography of the difference between conditions
                Odata(ic).Topography=TOP;
                [~,p,~,stats]=ttest(TOP(R_ind),TOP(L_ind));                                                             %Test Symmetry right vs left
                [top_value, top_index]=sort(abs((TOP-mean(TOP))./std(TOP)),'descend');                                  %Top 5 electrodes' distance from mean distribution?
            elseif ic==2
                TOP=squeeze(mean(mean(RegDat(:, tbin, I1),3),2))'-squeeze(mean(mean(RegDat(:, tbin, I2),3),2))';        %For BPall topo is topography of the difference between conditions
                Odata(ic).Topography=TOP;
                [~,p,~,stats]=ttest(TOP(R_ind),TOP(L_ind));                                                             %Test Symmetry right vs left
                [top_value, top_index]=sort(abs((TOP-mean(TOP))./std(TOP)),'descend');                                  %Top 5 electrodes' distance from mean distribution?
            else
                TOP=indata.icawinv(:,AC)';                                                                              %For ICs topo is topography of the component itself
                Odata(ic).Topography=TOP;
                [~,p,~,stats]=ttest(TOP(R_ind),TOP(L_ind));                                                             %Test Symmetry right vs left
                [top_value, top_index]=sort(abs((TOP-mean(TOP))./std(TOP)),'descend');                                  %Top 5 electrodes' distance from mean distribution?
            end;
            Odata(ic).top_value_abs         = TOP(1:5);
            Odata(ic).component_symmetry_p  = p;            %close to 1 the more symmetric
            Odata(ic).component_symmetry_t  = stats.tstat;  %t statistic of symmetry
            Odata(ic).top_lable             = {indata.chanlocs(top_index(1:5)).labels};
            Odata(ic).top_value_SD          = top_value(1:5);
            [~, ~, ~, ~, ~, outlier, outlier_num] = grubbs(TOP', {indata.chanlocs.labels});
            if ~isempty(outlier_num)
                OL=size(outlier_num,1); 
                if OL>5
                    OL=5;
                end;
                Odata(ic).outliers         = {outlier{1:OL,1}};
            else
                Odata(ic).outliers         = {'NaN'};
            end;

            %%%%%PERFORM REGRESSION IN FREQUENCY DOMAIN%%%%%%
            if ~isempty(Ffrequency)
                warning off
                    [Frequencies]=AGF_Hilbert(RegDat(AC, :, :), Ffreq_base, Ffrequency, [indata.xmin*1000 indata.xmax*1000], FTF_window(1), FTF_window(2), indata.srate,'divide');
                warning on
                tfrex=find(indata.times==FTF_window(1)):find(indata.times==FTF_window(2));
                interval=find(tfrex==Ftime_w(1)):find(tfrex==Ftime_w(end));
                for FQ = 1 : length(Ffrequency) %count through all specified frequency bands and save information to output structure
                    CurrFreq=Ffrequency{FQ};
                    ranks2 = nan(length(Ftime_w),12);                                
                    for bin = 1 : length(interval) %Rerun regression in frequency domain
                        %Define datapoints: Actual point + bin_size relative to time_window
                        tbin                    = [interval(bin)-ceil((Fbins-1)/2) : interval(bin)+floor((Fbins-1)/2)];
                        %Do this for all components
                        [t1, t2, ~, ~, STATS]   = regress(squeeze(mean(Frequencies(FQ,1,tbin,:),3)), [ones(length(trial_ind),1) trial_ind']);
                        ranks2(bin,1)           = STATS(2);    % remember F value
                        ranks2(bin,2)           = STATS(1);    % remember R2 value
                        ranks2(bin,3)           = STATS(3);    % remember p value
                        ranks2(bin,4)           = t1(2);       % remember B
                        ranks2(bin,5)           = t2(2,1);     % remember lower 95 CI
                        ranks2(bin,6)           = t2(2,2);     % remember upper 95 CI
                        ranks2(bin,7)           = mean(squeeze(mean(Frequencies(FQ,1,tbin,:),3)));                          % remember activity all
                        ranks2(bin,10)          = std(squeeze(mean(Frequencies(FQ,1,tbin,:),3)))/sqrt(length(trial_ind)-1); % remember std cat all
                        if length(vas)==2
                            ranks2(bin,8)       = mean(squeeze(mean(Frequencies(FQ,1,tbin,I1),3)));                         % remember activity cat 1
                            ranks2(bin,9)       = mean(squeeze(mean(Frequencies(FQ,1,tbin,I2),3)));                         % remember activity cat 2
                            ranks2(bin,11)      = std(squeeze(mean(Frequencies(FQ,1,tbin,I1),3)))/sqrt(length(I1)-1);       % remember std cat all
                            ranks2(bin,12)      = std(squeeze(mean(Frequencies(FQ,1,tbin,I2),3)))/sqrt(length(I2)-1);       % remember std cat all
                        end;
                    end;
                    Odata(ic).timecourse1_all     =ranks2(:,7);
                    Odata(ic).timecourse1_all_se  =ranks2(:,10);
                    if length(vas) == 2
                        Odata(ic).([ CurrLable  CurrFreq '_timecourse1'])=ranks2(:,8);
                        Odata(ic).([ CurrLable  CurrFreq '_timecourse2'])=ranks2(:,9);
                        Odata(ic).([ CurrLable  CurrFreq '_timecourse1_se'])=ranks2(:,11);
                        Odata(ic).([ CurrLable  CurrFreq '_timecourse2_se'])=ranks2(:,12);
                    else
                        Odata(ic).([ CurrLable  CurrFreq '_timecourse1'])=ranks2(:,8);
                        Odata(ic).([ CurrLable  CurrFreq '_timecourse2'])=NaN;
                        Odata(ic).([ CurrLable  CurrFreq '_timecourse1_se'])=ranks2(:,11);
                        Odata(ic).([ CurrLable  CurrFreq '_timecourse2_se'])=NaN;
                    end;
                    %Get stats to compare in the initial timewindow to select components
                    [x,y]=max(ranks2(Istat,1));y=y+Istat(1);
                    Odata(ic).([ CurrLable  CurrFreq '_Highest_F'])     = x;
                    Odata(ic).([ CurrLable  CurrFreq '_Highest_R2'])    = ranks2(y,2);
                    Odata(ic).([ CurrLable  CurrFreq '_Lowest_p'])      = ranks2(y,3);
                    Odata(ic).([ CurrLable  CurrFreq '_timepoint_ms'])  = indata.times([Ftime_w(1)+y])';
                    if isnan(Fbins)%Save best activity
                        tbin= Ftest_window+Ftime_w(1);
                    else
                        tbin = [interval(y)-ceil((Fbins-1)/2) : interval(y)+floor((Fbins-1)/2)];
                    end;
                    Odata(ic).([ CurrLable CurrFreq '_Value_Cat1'])            = mean(squeeze(mean(Frequencies(FQ,1,tbin,I1),3)));
                    Odata(ic).([ CurrLable CurrFreq '_Value_Cat1_se'])         = std(squeeze(mean(Frequencies(FQ,1,tbin,I1),3)))/sqrt(length(I1)-1);
                    Odata(ic).([ CurrLable CurrFreq '_Value_Cat2'])            = mean(squeeze(mean(Frequencies(FQ,1,tbin,I2),3)));
                    Odata(ic).([ CurrLable CurrFreq '_Value_Cat2_se'])         = std(squeeze(mean(Frequencies(FQ,1,tbin,I2),3)))/sqrt(length(I2)-1);
                    Odata(ic).([ CurrLable CurrFreq '_Activitiy_at_best_fit']) = squeeze(mean(Frequencies(FQ,1,tbin,:),3));
                    %save correlation coefficients and confidence intervals
                    Odata(ic).([ CurrLable CurrFreq '_beta_timecourse'])       = ranks2(:,4);
                    Odata(ic).([ CurrLable CurrFreq '_beta_timecourse_UCI'])   = ranks2(:,5);
                    Odata(ic).([ CurrLable CurrFreq '_beta_timecourse_LCI'])   = ranks2(:,6);
                end;
            end;
            %Get full spectral power in both conditions
            Fullfreq=[3 40];
            FreqSteps=39;
            Fspace='log';
            Fspectralbase=[-750 -500];
            Fnumcycles=4.5;
            if Fspeed~=1
                if ic == 1;disp(['Performing TF decomposition in ' num2str(Fullfreq(1)) ' to ' num2str(Fullfreq(2)) ' Hz, ' num2str(FreqSteps) ' steps, using ' Fspace '-space and ' num2str(Fnumcycles) ' cycles of Morlet wavelets.']);end;
                [SpectrumALL, ~, frex]=AGF_Tfreq( SpectDat, AC, 0, FreqSteps, Fullfreq(1), Fullfreq(2), Fspace, 4.5, FTF_window(1), FTF_window(2), Fspectralbase(1), Fspectralbase(2));
                SpectrumALL=SpectrumALL(:,:,2);
                Odata(ic).([ CurrLable 'Spectrum_ALL'])     = SpectrumALL;
                if length(vas)>1
                    UseData1=SpectDat;UseData1.data=SpectDat.data(:,:,I1);
                    [Spectrum1]=AGF_Tfreq( UseData1, AC, 0, FreqSteps, Fullfreq(1), Fullfreq(2), Fspace, 4.5, FTF_window(1), FTF_window(2), Fspectralbase(1), Fspectralbase(2) );
                    Spectrum1=Spectrum1(:,:,2);
                    UseData1=SpectDat;UseData1.data=SpectDat.data(:,:,I2);
                    [Spectrum2]=AGF_Tfreq( UseData1, AC, 0, FreqSteps, Fullfreq(1), Fullfreq(2), Fspace, 4.5, FTF_window(1), FTF_window(2), Fspectralbase(1), Fspectralbase(2) );
                    Spectrum2=Spectrum2(:,:,2);
                    Spectrum3=Spectrum1-Spectrum2;
                    Odata(ic).([ CurrLable 'Spectrum_1'])       = Spectrum1;
                    Odata(ic).([ CurrLable 'Spectrum_2'])       = Spectrum2;
                    Odata(ic).([ CurrLable 'Spectrum_DIFF'])    = Spectrum3;
                end;
            else
                Odata(ic).([ CurrLable 'Spectrum_1'])       = zeros(39,1001);
                Odata(ic).([ CurrLable 'Spectrum_2'])       = zeros(39,1001);
                Odata(ic).([ CurrLable 'Spectrum_DIFF'])    = zeros(39,1001);
            end;
            
            if nmodels>0 %model for more elaborate regression is provided!
                disp('Running single trial robust regression. Model no: ')
                for ModCount = 1 : nmodels
                    fprintf(2,'\b\b%d', ModCount);
                    UseModel = Fmodel{ModCount};
                    TI=UseModel{1}{:};
                    X =UseModel{2}{:};
                    Y =UseModel{3}{:};
                    Mname=UseModel{4}{:};

                    if size(X,2)>size(X,1)
                        X = X'; %transpose if dimensions do not match
                    end;
                    if size(Y,2)>size(Y,1)
                        Y = Y'; %transpose if dimensions do not match
                    end;

                    %Run that regression!
                    warning off
                    for Rcount = 1 : length(Ftime_w) 
                        tbin         = [Ftime_w(Rcount)-ceil((Fbins-1)/2) : Ftime_w(Rcount)+floor((Fbins-1)/2)];
                        [B, stats]   = robustfit([squeeze(mean(RegDat(AC, tbin, TI),2)) X], Y);
                        Tstat(Rcount)=stats.t(2);
                        Pstat(Rcount)=stats.p(2);
                    end;
                    warning on
                    Odata(ic).([ CurrLable 'Robust_Model' num2str(ModCount) '_Name'])=Mname;
                    Odata(ic).([ CurrLable 'Model_' num2str(ModCount) '_Timedomain_Robust_t' ]) = Tstat;
                    Odata(ic).([ CurrLable 'Model_' num2str(ModCount) '_Timedomain_Robust_p' ]) = Pstat;
                    Odata(ic).([ CurrLable 'Model_' num2str(ModCount) '_Timedomain_n_trials' ]) = length(TI);
                    [~,y]=max(abs(Tstat(Istat)));
                    y=y+Istat(1)-1;
                    Odata(ic).([ CurrLable 'Model_' num2str(ModCount) '_Timedomain_Highest_t' ]) = Tstat(y);
                    Odata(ic).([ CurrLable 'Model_' num2str(ModCount) '_Timedomain_Lowest_p' ])  = Pstat(y);
                    Odata(ic).([ CurrLable 'Model_' num2str(ModCount) '_Timedomain_timepoint_ms' ])  = indata.times([Ftime_w(1)+y])';


                    %Run that regression in frequency domain (only the first frequency is used)!
                    if ~isempty(Ffrequency)
                        CurrFreq=Ffrequency{1};
                        warning off
                        for Rcount = 1 : length(interval)
                            tbin                    = [interval(Rcount)-ceil((Fbins-1)/2) : interval(Rcount)+floor((Fbins-1)/2)];
                            [B, stats]   = robustfit([squeeze(mean(RegDat(AC, tbin, TI),2)) X], Y);
                            Tstat(Rcount)=stats.t(2);
                            Pstat(Rcount)=stats.p(2);
                        end;
                        warning on 
                        Odata(ic).([ CurrLable 'Model_' num2str(ModCount) '_' CurrFreq '_Robust_t' ]) = Tstat;
                        Odata(ic).([ CurrLable 'Model_' num2str(ModCount) '_' CurrFreq '_Robust_p' ]) = Pstat;
                        Odata(ic).([ CurrLable 'Model_' num2str(ModCount) '_' CurrFreq '_n_trials' ]) = length(TI);
                        [~,y]=max(abs(Tstat(Istat)));
                        y=y+Istat(1)-1;
                        Odata(ic).([ CurrLable 'Model_' num2str(ModCount) '_' CurrFreq '_Highest_t' ]) = Tstat(y);
                        Odata(ic).([ CurrLable 'Model_' num2str(ModCount) '_' CurrFreq '_Lowest_p' ]) = Pstat(y);
                        Odata(ic).([ CurrLable 'Model_' num2str(ModCount) '_' CurrFreq '_timepoint_ms' ])  = indata.times([Ftime_w(1)+y])';
                    end;
                fprintf(2,'\n')
                end;
            end;
        end;
    end;
else
    Odata=[];
end;
add_info.icawinv=indata.icawinv;
add_info.icasphere=indata.icasphere;
add_info.chanlocs=indata.chanlocs;
return;



































